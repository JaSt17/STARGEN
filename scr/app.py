import streamlit as st
import folium
from streamlit_folium import folium_static
import pandas as pd
import os
from vizualize import *
from func import *

def clear_state():
    """
    Clear all keys from the Streamlit session state.
    """
    for key in list(st.session_state.keys()):
        del st.session_state[key]
        

def get_resolution_data():
    """
    Create and return a DataFrame containing resolution data for display.
    
    Returns:
        pd.DataFrame: DataFrame with resolution information.
    """
    data = {
        "Resolution": [0, 1, 2, 3, 4],
        "Total number of cells": [122, 842, 5882, 41162, 288122],
        "Average cell area (km2)": [4357449.41, 609788.44, 86801.78, 12393.43, 1770.34],
        "Average edge length (Km)": [1281.256011, 483.0568391, 182.5129565, 68.97922179, 26.07175968]
    }
    df = pd.DataFrame(data)
    return df


def initialize_session():
    """
    Initialize the session state and display the initial setup UI.
    """
    st.set_page_config(page_title="STARGEN (Spatio-Temporal Analysis and Reconstruction of GENetic Barriers)", page_icon=":earth_americas:")
    
    # Display the logo and title
    col1, col2 = st.columns([1, 2])
    # change the columns width
    with col1:
        st.title('STARGEN')
    with col2:
        st.image("img/STARGEN.png", width=80)
    st.write('**(Spatio-Temporal Analysis and Reconstruction of GENetic Barriers)**')

    # Slider for selecting the number of time bins
    st.session_state['time_bins'] = st.slider('Select a number of time bins', 5, 30, 14, 1)
    # Checkbox to enable same time bin length
    st.session_state['same_age_range'] = st.checkbox('Same age range for each time bin', value=True)

    # Display the current time span for each time bin
    if st.session_state['same_age_range']:
        st.write(f"Current time span for each time bin is about {round((14000)/10/st.session_state['time_bins'])*10} years.")

    # Button to display information about time bins
    if st.button("Information about time bins"):
        st.write("""
        The ancient samples' dates range from 1890 AD to 108500 BC.
        The data will get organized into time bins, which hold the same number of samples per time bin.
        The samples can also be arranged so every time bin has the same range of years by checking the checkbox above.
        In that case the last timebin will be longer than the others to include all samples.
        """)

    # Slider for selecting the resolution
    st.session_state['resolution'] = st.slider('Select a resolution', 1, 4, 3, 1)

    # Button to display information about resolution
    if st.button("Information about resolution"):
        st.table(get_resolution_data())

    # Button to run the tool
    if st.button('Run'):
        st.text('Running STARGEN...')
        st.session_state['setup_done'] = True
        # Load the distance matrix to the session state
        path_to_matrix = os.getcwd() + "/1_dist_matrix/eucl_dist.pkl"
        st.session_state['matrix'] = pd.read_pickle(path_to_matrix)
        st.rerun()
        

def setup_done_ui():
    """
    Display the main UI after the setup is done.
    """
    # Return to Home button
    if st.button('Return to Home', key='home'):
        clear_state()
        st.rerun()

    # Label the samples and save them in a DataFrame
    if 'df' not in st.session_state:
        st.session_state['df'] = label_samples(os.getcwd(), st.session_state['time_bins'], st.session_state['resolution'], st.session_state['same_age_range'])
    
    # Calculate the average distances between neighboring hexagons for each time bin
    if 'time_bins_dist' not in st.session_state:
        st.session_state['time_bins_dist'] = calc_dist_time_bin(st.session_state['df'], st.session_state['matrix'])
    
    # Rename the time bins to display them in the dropdown
    time_bins = rename_time_bins(st.session_state['df'])
    selected_time_bin = st.selectbox("Time Bin", options=time_bins)

    if 'selected_time_bin_id' not in st.session_state:
        st.session_state['selected_time_bin_id'] = 1
    
    # get the id of the selected time bin
    new_selected_time_bin_id = time_bins.index(selected_time_bin)
    # get the time bin and the hexagons for the selected time bin
    time_bin = st.session_state['time_bins_dist'][selected_time_bin]
    # scale the distances to their geographical distance and save the predicted distances to scale the internal distances with it
    time_bin, gen_distances_pred = scale_distances(time_bin, resolution=st.session_state['resolution'])
    # get the hexagons with there internal distance and the distance values for the selected time bin
    time_bin, hexagons = get_hexagons(time_bin)

    # Initialize thresholds for distance values and isolated populations
    if 'threshold' not in st.session_state:
        st.session_state['threshold'] = -5.0
        st.session_state['isolated_threshold'] = 1.0

    # Slider to choose the threshold for the distance values to display
    st.session_state['threshold'] = st.sidebar.slider('Minimal distance value to display?', -5.0, 5.0, st.session_state['threshold'], 0.1)
    
    # Slider to choose the threshold for isolated populations
    new_isolated_threshold = st.sidebar.slider('Minimal distance value to be considered as isolated?', 0.0, 4.0, st.session_state['isolated_threshold'], 0.1)
    
    new_allowed_distance = st.sidebar.slider('Number of hexagons to consider as neighborhood?', 1, 25, 10)

    # Check if the threshold or the selected time bin has changed
    if new_isolated_threshold != st.session_state['isolated_threshold'] or new_selected_time_bin_id != st.session_state['selected_time_bin_id'] or new_allowed_distance != st.session_state['allowed_distance']:
        # save the new allowed distance 
        st.session_state['allowed_distance'] = new_allowed_distance
        # get the isolated hexagons and the barrier lines for the selected time bin
        st.session_state['isolated_hex'], st.session_state['barrier_lines'], st.session_state['barrier_hex'], st.session_state['new_time_bin'] = get_isolated_hex_and_barriers(time_bin, hexagons, st.session_state['isolated_threshold'], st.session_state['allowed_distance'])
        # impute the missing hexagons until a range of 2 times the resolution is reached
        st.session_state['imputed_hex'] = impute_missing_hexagons(st.session_state['barrier_hex'], num_runs=st.session_state['resolution'] * 2)
        # save the new threshold for isolated populations
        st.session_state['isolated_threshold'] = new_isolated_threshold
        # save the new selected time bin id
        st.session_state['selected_time_bin_id'] = new_selected_time_bin_id
        # find the closest populations to the isolated populations
        st.session_state['closest_populations'], st.session_state['isolated_hex'] = find_closest_population(
            st.session_state['df'], st.session_state['selected_time_bin_id'], st.session_state['isolated_hex'], 
            st.session_state['matrix'], st.session_state['isolated_threshold'], gen_distances_pred, st.session_state['resolution'])

    # Checkbox to toggle showing possible migration routes
    st.session_state['show_migration'] = st.sidebar.checkbox("Show possible migration routes", False)
    # Checkbox to toggle showing isolated populations
    st.session_state['show_isolated'] = st.sidebar.checkbox("Show isolated populations", False)
    # Checkbox to toggle showing distance lines
    st.session_state['show_lines'] = st.sidebar.checkbox("Show line representation", False)

    # Set initial map state if it does not exist
    if 'map_state' not in st.session_state:
        st.session_state['map_state'] = {"lat": 42.0, "lon": 44.75, "zoom": 1}

    st.sidebar.write("Specify a default window for the map:")
    st.session_state['map_state']['lat'] = st.sidebar.number_input('Enter latitude:', -90.0, 90.0, step=0.01, value=st.session_state['map_state']['lat'])
    st.session_state['map_state']['lon'] = st.sidebar.number_input('Enter longitude:', -180.0, 180.0, step=0.01, value=st.session_state['map_state']['lon'])
    st.session_state['map_state']['zoom'] = st.sidebar.slider('Choose zoom level:', 1, 10, st.session_state['map_state']['zoom'])

    lat, lon, zoom = st.session_state['map_state'].values()
    map_tiles = "Esri worldstreetmap" if not st.sidebar.checkbox("Black and White Map", False) else "Cartodb Positron"
    m = folium.Map(location=(lat, lon), tiles=map_tiles, zoom_start=zoom)

    # Draw lines or hexagons based on the selected options
    if st.session_state['show_lines']:
        m = draw_sample_hexagons(hexagons, m, zoom_start=zoom)
        # use the new time bin to only show distnaces that are in the allowed distance
        lines = get_distance_lines(st.session_state['new_time_bin'])
        m = draw_barriers(lines, m, threshold=-10.0)
    else:
        m = draw_hexagons_with_values(st.session_state['barrier_hex'], m, threshold=st.session_state['threshold'])
        m = draw_hexagons_with_values(st.session_state['imputed_hex'], m, threshold=st.session_state['threshold'], imputed=True)
        m = draw_sample_hexagons(hexagons, m, zoom_start=zoom)

    # Draw migration routes if selected
    if st.session_state['show_migration']:
        m = draw_migration_for_time_bin(st.session_state['closest_populations'], m)
        m = draw_sample_hexagons(hexagons, m, zoom_start=zoom)
    
    # Draw isolated populations if selected
    if st.session_state['show_isolated']:
        m = draw_hexagons(st.session_state['isolated_hex'], m, color="black", opacity=0.6)
        m = draw_sample_hexagons(hexagons, m, zoom_start=zoom)
        
    # Draw barriers if there are any
    if len(st.session_state['barrier_lines']) > 0:
        m = draw_barriers(st.session_state['barrier_lines'], m, threshold=st.session_state['threshold'])
        
    m = add_legend(m)
    folium_static(m, width=800, height=600)
    st.write(f"Number of isolated populations: {len(st.session_state['isolated_hex'])}")
    

def main():
    """
    Main function to control the app's flow.
    """
    if 'setup_done' not in st.session_state:
        initialize_session()
    else:
        setup_done_ui()

if __name__ == "__main__":
    main()
