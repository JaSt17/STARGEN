import streamlit as st
import folium
from streamlit_folium import folium_static
import pandas as pd
import pickle
import os
from label_samples_time_hexa import label_samples
from vizualize import draw_hexagons, draw_all_boarders_for_time_bin, add_colorbar_to_map
from func import rename_time_bins, calc_dist_time_bin, normalize_distances, get_time_bin_hexagons, get_min_max_dist

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
    time_bins = st.slider('Select a number of time bins', 1, 30, 18, 1)
    # Checkbox to enable same time bin length
    same_age_range = st.checkbox('Same age range for each time bin', value=True)

    # show the user the current time span for each time bin give the number of time bins
    if same_age_range:
        st.write(f"Current time span for each time bin is {round((18000)/time_bins)} years.")
        
    # button to get information about time bins
    if st.button("Information about time bins"):
            st.write("""
    The ancient samples' dates range from 1890 AD to 108500 BC.
    The data will get organized into time bins, which hold the same number of samples per time bin.
    The samples can also be arrange so every time bin has the same range of years by checking the checkbox above.
        """)
    
    # Slider to choose the resolution and button to get information about it
    resolution = st.slider('Select a resolution', 0, 5, 3, 1)
    
    # button to get information about resolution
    if st.button("Information about resolution"):
        st.table(get_resolution_data())

    # Checkbox to allow expanding the search area and button to get information about it
    allow_k_distance = st.checkbox('Expand Search Area', value=True)
    
    # button to get information about expanding search area
    if st.button("Information about search area & neighbors"):
        st.write("""
    As the resolution increases, the size of each hexagon decreases, resulting in a reduced likelihood of adjacent hexagons being direct neighbors. 
    To ensure visibility between closely situated hexagons, users are encouraged to utilize the "Expand Search Area" option.
    When activated, this feature systematically searches through hexagons at greater distances until a neighboring hexagon is identified.
    """)
    
    # Button to run the tool
    if st.button('Run'):
        st.text('Running GeoGeneTrack...')
        # Save the selected parameters to the session state
        st.session_state['setup_done'] = True
        st.session_state['same_age_range'] = same_age_range
        st.session_state['time_bins'] = time_bins
        st.session_state['resolution'] = resolution
        st.session_state['same_age_range'] = st.session_state['same_age_range']
        st.session_state['allow_k_distance'] = allow_k_distance
        # read the distance matrix from the file
        # get path to the distance matrix
        path_to_matrix = os.getcwd()+"/1_dist_matrix/eucl_dist.pkl"
        st.session_state['matrix'] = pd.read_pickle(path_to_matrix)
        st.rerun()
        
# Once setup is done, show the map and time bin selection 
if 'setup_done' in st.session_state and st.session_state['setup_done']:
    # Return to Home button
    if st.sidebar.button('Return to Home', key='home'):
        clear_state()
        st.rerun()
        
    # initialize the neighborhood size for the distance calculation and the scale by distance
    if 'neighborhood_size' not in st.session_state:
        st.session_state['neighborhood_size'] = 1
        st.session_state['scale_by_distance'] = False
        
    @st.cache_data
    def load_data():
        return label_samples(os.getcwd(),st.session_state['time_bins'],st.session_state['resolution'], st.session_state['same_age_range'])
    df = load_data()
    
    @st.cache_data
    def load_time_bins():
        return rename_time_bins(df)
    time_bins = load_time_bins()
    
    @st.cache_data
    def load_hexagons():
        return get_time_bin_hexagons(df)
    time_bins_hexagons = load_hexagons()

    # dropdown to select time bins
    selected_time_bin = st.sidebar.selectbox("Time Bin", options=time_bins)
    
    # checkbox for normalizing distance values
    if st.sidebar.checkbox("Normalize distances", value=False):
        # normalize the distance values
        time_bin = normalize_distances(time_bin)
        
    # checkbox for scaling the distance values by the distance
    if st.sidebar.checkbox("Scale distances by distance", value=False):
        st.session_state['scale_by_distance'] = True
    else:
        st.session_state['scale_by_distance'] = False
        
    # calculate the average distances between neighboring hexagons for each time bin with the given parameters
    time_bins_dist = calc_dist_time_bin(df, st.session_state['matrix'],
                                        st.session_state['neighborhood_size'] ,
                                        st.session_state['allow_k_distance'],
                                        st.session_state['scale_by_distance'])

    # get the hexagons and distance values for the selected time bin
    hexagons = time_bins_hexagons[selected_time_bin]
    time_bin = time_bins_dist[selected_time_bin]
    
    # get min and max distance for the selected time bin for the threshold
    min_dist, max_dist = get_min_max_dist(time_bin)
    
    # initialize the threshold with the minimum distance value
    if 'threshold' not in st.session_state:
        st.session_state['threshold'] = min_dist
        
    # change the threshold if the user changes the time bin
    if st.session_state['threshold'] < min_dist:
        st.session_state['threshold'] = min_dist
    if st.session_state['threshold'] > max_dist:
        st.session_state['threshold'] = max_dist
        
    # Slider to choose the threshold for the distant values which should be displayed
    st.session_state['threshold'] = st.sidebar.slider('choose threshold:', min_dist, max_dist, st.session_state['threshold'], 0.01)
    
    # set an initial map state if it does not exist
    if 'map_state' not in st.session_state:
        st.session_state['map_state'] = {
            "lat": 42.0,
            "lon": 44.75,
            "zoom": 2
        }
    st.sidebar.write("You can change the default latitude, longitude, and zoom level.")
    # text input where the user can enter the latitude between -90 and 90
    st.session_state['map_state']['lat'] = st.sidebar.number_input('Enter latitude:', -90.0, 90.0, step=0.01, value=st.session_state['map_state']['lat'])
    # text input where the user can enter the longitude between -180 and 180
    st.session_state['map_state']['lon'] = st.sidebar.number_input('Enter longitude:', -180.0, 180.0, step=0.01, value=st.session_state['map_state']['lon'])
    # slider to choose the zoom level of the map
    st.session_state['map_state']['zoom'] = st.sidebar.slider('Choose zoom level:', 1, 15, st.session_state['map_state']['zoom'])
    
    # initialize the neighborhood size for the distance calculation
    if 'neighborhood_size' not in st.session_state:
        st.session_state['neighborhood_size'] = 1
        
    # Slider to choose the neighborhood size for the distance calculation
    st.session_state['neighborhood_size'] = st.sidebar.slider('Choose neighborhood size:', 1, 10, st.session_state['neighborhood_size'])
    
    if button := st.sidebar.button('What is the neighborhood size?'):
        st.write("""The number of neighbors to search can be adjusted using the slider above. The default value is 1. If we increase the number of neighbors to search,
    the tool will search for more neighbors until reacing the selected neighborhood instead of stoping when the first neighbor is found.""")

    if 'map_state' in st.session_state:
        lat, lon, zoom = st.session_state['map_state'].values()
        m = folium.Map(location=(lat, lon), tiles="Esri worldstreetmap", zoom_start=zoom)
        m = draw_hexagons(hexagons, m, zoom_start=zoom)
    else:
        m = draw_hexagons(hexagons)
        
    # draw all average distance values between neighboring hexagons for the chosen time bin
    m = draw_all_boarders_for_time_bin(time_bin, m, threshold=st.session_state['threshold'])
    # Display the current time bin and chosen threshold on top the map
    st.markdown(f"<h1 style='text-align: center; font-size: 20px;'>Time Bin: {selected_time_bin}</h1>", unsafe_allow_html=True)
    # Display the map in Streamlit
    folium_static(m, width=800, height=600)