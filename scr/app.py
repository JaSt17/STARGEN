import streamlit as st
import folium
from streamlit_folium import folium_static
import pandas as pd
import pickle
import os
from label_samples_time_hexa import label_samples
from vizualize import draw_hexagons, draw_migration_for_time_bin, draw_hexagons_with_values, draw_barriers
from func import *

# Function to clear session state
def clear_state():
    for key in list(st.session_state.keys()):
        del st.session_state[key]
            
# Function to get resolution data as table to display it for the user
def get_resolution_data():
    data = {
    "Resolution": [0, 1, 2, 3, 4, 5],
    "Total number of cells": [122, 842, 5882, 41162, 288122, 2016842],
    "Average cell area (km2)": [4357449.41, 609788.44, 86801.78, 12393.43, 1770.34, 252.90],
    "Average edge length (Km)":  [1281.256011, 483.0568391, 182.5129565, 68.97922179, 26.07175968, 9.854090990]
    }
    # Convert the dictionary to a pandas DataFrame
    df = pd.DataFrame(data)
    return df

# Initial screen to select time bins and resolution
if 'setup_done' not in st.session_state:
    st.set_page_config(page_title="GeoGeneTrack", page_icon=":earth_americas:")
    # display the logo and title
    col1, col2 = st.columns([1, 2])
    with col2:
        st.image("img/GeoGenTrack_logo.png", width=100) 
    with col1:
        st.title('GeoGenTrack')
    
    # Slider to choose the number of time bins and button to get information about it
    st.session_state['time_bins'] = st.slider('Select a number of time bins', 1, 20, 11, 1)
    # Checkbox to enable same time bin length
    st.session_state['same_age_range'] = st.checkbox('Same age range for each time bin', value=True)

    # show the user the current time span for each time bin give the number of time bins
    if st.session_state['same_age_range']:
        st.write(f"Current time span for each time bin is {round((11000)/st.session_state['time_bins'])} years.")
        
    # button to get information about time bins
    if st.button("Information about time bins"):
            st.write("""
    The ancient samples' dates range from 1890 AD to 108500 BC.
    The data will get organized into time bins, which hold the same number of samples per time bin.
    The samples can also be arrange so every time bin has the same range of years by checking the checkbox above.
        """)
    
    # Slider to choose the resolution and button to get information about it
    st.session_state['resolution'] = st.slider('Select a resolution', 0, 5, 3, 1)
    
    # button to get information about resolution
    if st.button("Information about resolution"):
        st.table(get_resolution_data())
        
    # Slider to choose the neighborhood size for the distance calculation
    st.session_state['neighborhood_size'] = st.slider('Neighborhood range:', 1, 20, 7)
    
        # button to get information about time bins
    if st.button("Information about Neighborhood range"):
            st.write("""
    The neighborhood range is the number of hexagons that are considered as neighbors for each hexagon. For example, if the neighborhood range is 5,
    the average genetic distance between a hexagon and all hexagons in a range of 5 hexagons will be calculated. Increasing the neighborhood range will increase the computation time.
        """)

    # Button to run the tool
    if st.button('Run'):
        st.text('Running GeoGeneTrack...')
        # Save the selected parameters to the session state
        st.session_state['setup_done'] = True
        # read the distance matrix from the file
        # get path to the distance matrix
        path_to_matrix = os.getcwd()+"/1_dist_matrix/eucl_dist.pkl"
        # load the distance matrix to the session state
        st.session_state['matrix'] = pd.read_pickle(path_to_matrix)
        st.rerun()
        
# Once setup is done, show the map and time bin selection 
if 'setup_done' in st.session_state and st.session_state['setup_done']:
    # Return to Home button
    if st.button('Return to Home', key='home'):
        clear_state()
        st.rerun()
    
    if 'df' not in st.session_state:
        # lable the samples with hexagon id and time bin and save it in a dataframe
        st.session_state['df'] = label_samples(os.getcwd(),st.session_state['time_bins'],st.session_state['resolution'], st.session_state['same_age_range'])
    if 'time_bins_dist' not in st.session_state:
        # calculate the average distances between neighboring hexagons for each time bin with the given parameters
        st.session_state['time_bins_dist'] = calc_dist_time_bin(st.session_state['df'], st.session_state['matrix'], st.session_state['neighborhood_size'], False)
    
    # rename the time bins to display them in the dropdown
    time_bins = rename_time_bins(st.session_state['df'])
    # get the hexagons for each time bin
    time_bins_hexagons = get_time_bin_hexagons(st.session_state['df'])

    # dropdown to select time bins
    selected_time_bin = st.selectbox("Time Bin", options=time_bins)
    if 'selected_time_bin_id' not in st.session_state:
        st.session_state['selected_time_bin_id'] = 1
    new_selected_time_bin_id = time_bins.index(selected_time_bin)

    # get the hexagons and distance values for the selected time bin
    hexagons = time_bins_hexagons[selected_time_bin]
    # get the distance values for the selected time bin
    time_bin = st.session_state['time_bins_dist'][selected_time_bin]
    # normalize the distances for each timebin
    time_bin = normalize_distances(time_bin)
    
    # initialize the threshold for the distance values to display and the threshold for isolated populations
    if 'threshold' not in st.session_state:
        st.session_state['threshold'] = 0.0
        st.session_state['isolated_threshold'] = 0.4
        
    # Slider to choose the threshold for the distant values which should be displayed
    st.session_state['threshold'] = st.sidebar.slider('Which distances should be displayed ?:',
                                                        0.0, 1.0, st.session_state['threshold'], 0.01)
    
    # Slider to chose the threshold for isolated populations
    new_isolated_threshold = st.sidebar.slider('Which distances are considered as isolated populations ?:',
                                                0.0, 1.0, st.session_state['isolated_threshold'], 0.01)
    
    # check if the threshold or the selected time bin has changed
    if new_isolated_threshold != st.session_state['isolated_threshold'] or new_selected_time_bin_id != st.session_state['selected_time_bin_id']:
        # get the isolated hexagons and barriers for the selected time bin
        st.session_state['isolated_hex'], st.session_state['barrier_lines'], st.session_state['barrier_hex'] = get_isolated_hex_and_barriers(time_bin, hexagons, st.session_state['isolated_threshold'])
        # get imputed hexagons
        st.session_state['imputed_hex'] = impute_missing_hexagons_multiple_runs(st.session_state['barrier_hex'], hexagons, num_runs=8)
        # change the threshold for isolated populations
        st.session_state['isolated_threshold'] = new_isolated_threshold
        # change the selected time bin id
        st.session_state['selected_time_bin_id'] = new_selected_time_bin_id
        # find the closest populations for the isolated hexagons
        st.session_state['closest_populations'], st.session_state['isolated_hex'] = find_closest_population(st.session_state['df'],
                                                                                                            st.session_state['selected_time_bin_id'],
                                                                                                            st.session_state['isolated_hex'],
                                                                                                            st.session_state['matrix'],
                                                                                                            st.session_state['isolated_threshold'])
    
    if st.sidebar.checkbox("Show possible migration routes", False):
        st.session_state['show_migration'] = True
    else:
        st.session_state['show_migration'] = False
    
    # set an initial map state if it does not exist
    if 'map_state' not in st.session_state:
        st.session_state['map_state'] = {
            "lat": 42.0,
            "lon": 44.75,
            "zoom": 1
        }
    st.sidebar.write("Specify a default window for the map:")
    # text input where the user can enter the latitude between -90 and 90
    st.session_state['map_state']['lat'] = st.sidebar.number_input('Enter latitude:', -90.0, 90.0, step=0.01, value=st.session_state['map_state']['lat'])
    # text input where the user can enter the longitude between -180 and 180
    st.session_state['map_state']['lon'] = st.sidebar.number_input('Enter longitude:', -180.0, 180.0, step=0.01, value=st.session_state['map_state']['lon'])
    # slider to choose the zoom level of the map
    st.session_state['map_state']['zoom'] = st.sidebar.slider('Choose zoom level:', 1, 15, st.session_state['map_state']['zoom'])
        
    # get the specifyed map state
    lat, lon, zoom = st.session_state['map_state'].values()
    if not st.sidebar.checkbox("Black and White Map", False):
        m = folium.Map(location=(lat, lon), tiles="Esri worldstreetmap", zoom_start=zoom)
    else:
        m = folium.Map(location=(lat, lon),  tiles="Cartodb Positron", zoom_start=zoom)
    # draw all hexagons for the selected time bin which hold the samples
    m = draw_hexagons(hexagons, m, zoom_start=zoom)
    # check if the migration routes should be displayed
    if st.session_state['show_migration']:
        # draw the migration lines for the isolated hexagons
        m = draw_migration_for_time_bin(st.session_state['closest_populations'], m)
    # highlight the isolated hexagons in red that can not be explained by migration
    m = draw_hexagons(st.session_state['isolated_hex'], m, color='red', zoom_start=zoom, opacity=0.7, value='Isolated Population without migration route.')
    # draw the hexagons barriers and barrier lines between direct neighbors
    m = draw_hexagons_with_values(st.session_state['barrier_hex'], m, threshold = st.session_state['threshold'])
    # draw the imputed hexagons
    m = draw_hexagons_with_values(st.session_state['imputed_hex'], m, threshold = st.session_state['threshold'], imputed=True)
    # check if there are any barriers
    if len(st.session_state['barrier_lines']) > 0:
        m = draw_barriers(st.session_state['barrier_lines'], m, threshold = st.session_state['threshold'])
    # Display the map in Streamlit
    folium_static(m, width=800, height=600)
    st.write(f"Number of isolated populations with no migration route found: {len(st.session_state['isolated_hex'])}")