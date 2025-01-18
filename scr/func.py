import pandas as pd
import math
import numpy as np
import h3
from collections import defaultdict
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
import statsmodels.api as sm
from haversine import haversine

def read_df(path):
    """
    Reads the ancient DNA Annotations file into a DataFrame.

    Args:
        path (str): Path to the ancient DNA annotations file.

    Returns:
        pd.DataFrame: DataFrame containing the ancient DNA annotations with the 'Age' column converted to numeric.
    
    """
    # Read the file into a DataFrame
    df = pd.read_csv(path, sep="\t")
    # Check if the 'Age' column exists in the DataFrame
    if 'Age' not in df.columns:
        raise KeyError(f"The 'Age' column is missing from the file at {path}.")
    # Convert the 'Age' column to numeric, coercing errors to NaN
    df['Age'] = pd.to_numeric(df['Age'], errors='coerce')
    
    return df


def calc_avg_dist(samples_hex1, samples_hex2, dist_matrix):
    """
    Calculates the average distance between two groups of samples.

    Args:
        samples_hex1 (list): List of sample identifiers for the first group.
        samples_hex2 (list): List of sample identifiers for the second group.
        dist_matrix (pd.DataFrame): DataFrame containing the distance matrix.
    
    Returns:
        float: The average distance between the two groups of samples. 
    """
    return dist_matrix.loc[samples_hex1, samples_hex2].values.flatten().mean()


def calc_neighbor_dist(hexagons, dist_matrix, time_bin_df, hex_col):
    """
    Calculate the average distances between neighboring hexagons.

    Parameters:
    - hexagons (list): List of hexagon IDs.
    - dist_matrix (pd.DataFrame): A DataFrame representing the distance matrix.
    - time_bin_df (pd.DataFrame): A DataFrame containing the Data for the time_bin with a column for hexagon IDs.
    - hex_col (str): Column name in time_bin_df that contains hexagon IDs.

    Returns:
    - dict: A dictionary where keys are pairs of hexagons (as frozensets) and values are the average distances between them.
    """

    # Get the samples in each hexagon
    samples_in_hex = time_bin_df.groupby(hex_col)['ID'].apply(list).to_dict()
    
    # Create a list of all samples in the hexagons
    all_samples = [sample for samples in samples_in_hex.values() for sample in samples]
    
    # Create a submatrix of the distance matrix for the samples in the hexagons
    dist_matrix = dist_matrix.loc[all_samples, all_samples]
    
    # Initialize the dictionary to store the average distances between neighboring hexagons
    averages = {}

    def get_neighbors(hexagons):
        """
        Get neighbors for all hexagons using Delaunay triangulation.

        Parameters:
        - hexagons (list): List of hexagon IDs.

        Returns:
        - dict: A dictionary where keys are hexagon IDs and values are lists of neighboring hexagon IDs.
        """
        # Get the centroid of each hexagon
        coords = np.array([h3.cell_to_latlng(hex) for hex in hexagons])
        hexagons = np.array(hexagons)
        
        # Calculate the Delaunay triangulation
        tri = Delaunay(coords)

        # Get the edges (pairs of hexagons)
        edges = set()
        for simplex in tri.simplices:
            for i in range(3):
                edge = tuple(sorted((simplex[i], simplex[(i + 1) % 3])))
                edges.add(edge)
        
        # Create the dictionary of neighbors
        neighbors = {}
        for i, j in edges:
            neighbors.setdefault(hexagons[i], []).append(hexagons[j])
            neighbors.setdefault(hexagons[j], []).append(hexagons[i])
        
        return neighbors

    # Get the neighbors for whom we will calculate the distances
    neighbors = get_neighbors(hexagons)

    # Add hexagons not yet in the neighbors dictionary with an empty list
    for hexagon in hexagons:
        if hexagon not in neighbors:
            neighbors[hexagon] = []

    # Calculate the average distance between the hexagon and its neighbors
    for hexagon, neighbor_list in neighbors.items():
        # append the hexagon to the neighbor list to calculate the distance with itself
        neighbor_list.append(hexagon)
        for neighbor in neighbor_list:
            ids_in_hexagon = samples_in_hex.get(hexagon, [])
            ids_in_neighbor = samples_in_hex.get(neighbor, [])
            # Get the pair of hexagons
            pair = frozenset([hexagon, neighbor])

            # Calculate the average distance between the hexagon and its neighbor
            distance = calc_avg_dist(ids_in_hexagon, ids_in_neighbor, dist_matrix)

            averages[pair] = round(distance, 5)

    return averages


def calc_dist_time_bin(df, dist_matrix=None):
    """
    Calculate the average distance between each hexagon and its neighbors for each time bin.

    Parameters:
    - df (pd.DataFrame): DataFrame containing the data with hexagon IDs and age groups.
    - dist_matrix (pd.DataFrame): A DataFrame representing the distance matrix. Default is None.

    Returns:
    - dict: A dictionary where keys are time bin labels and values are dictionaries of average distances 
            between neighboring hexagons for each time bin.
    """
    
    # Get the column name for hexagons (it should be the only column with 'hex' in the name)
    hex_col = df.columns[df.columns.str.contains('hex')][0]
    
    # Convert 'AgeGroup' column values to tuples of integers representing the start and end years
    df['AgeGroupTuple'] = df['AgeGroup'].apply(lambda x: tuple(map(int, x.split('-'))))
    
    # Sort the unique age group tuples to process them in chronological order
    time_bins = sorted(df['AgeGroupTuple'].unique())
    averages = {}
    number_of_samples = {}

    # Iterate over each time bin
    for time_bin in time_bins:
        # Format the current time bin as a string for labeling purposes
        bin_label = rename_times(time_bin)
        
        # Get subset of the DataFrame for the current time bin
        time_bin_df = df[df['AgeGroupTuple'] == time_bin]
        
        # Get the number of samples in each hexagon
        samples_in_hex = time_bin_df.groupby(hex_col)['ID'].apply(list).to_dict()
        samples_in_hex = {hex: len(samples) for hex, samples in samples_in_hex.items()}
        
        # Append the number of samples in each hexagon to the dictionary using the time bin label as the key
        number_of_samples[bin_label] = samples_in_hex

        # Get all unique hexagons for the current time bin
        hexagons = time_bin_df[hex_col].unique()
        
        # Calculate the average distance for each hexagon to its neighbors within the current time bin
        average_distances = calc_neighbor_dist(hexagons, dist_matrix, time_bin_df, hex_col)

        # Append the calculated average distances to the dictionary using the time bin label as the key
        averages[bin_label] = average_distances
        

    # Return the dictionary with the average distances between neighboring hexagons for each time bin
    return averages, number_of_samples


def get_hexagons(time_bin):
    """
    Separates the hexagons from a time bin dictionary.

    Args:
        time_bin (dict): Dictionary where keys are pairs of samples and values are internal distances.

    Returns:
        tuple: A tuple containing:
            - new_time_bin (dict): Dictionary with pairs of samples as keys.
            - hexagons (dict): Dictionary with single samples as keys.
    """
    # Separate entries with single samples from those with pairs of samples
    hexagons = {list(pair)[0]: distance for pair, distance in time_bin.items() if len(pair) == 1}
    new_time_bin = {pair: distance for pair, distance in time_bin.items() if len(pair) > 1}
    
    return new_time_bin, hexagons


def get_isolated_hex_and_barriers(time_bin, hexagons, threshold, allowed_distance=12):
    """
    Identify isolated hexagons and barriers for each time bin.

    Parameters:
    - time_bin (dict): Dictionary containing pairs of hexagons and their distances.
    - hexagons (list): List of hexagon IDs.
    - threshold (float): Distance threshold to consider a hexagon isolated.
    - allowed_distance (int, optional): Maximum allowed distance for a line of hexagons. Default is 12.

    Returns:
    - tuple: 
        - isolated_hex (list): List of hexagons considered isolated.
        - barrier_lines (dict): Dictionary of shared boundaries with their distances.
        - barrier_hex (dict): Dictionary of hexagons with their average barrier distances.
        - new_time_bin (dict): new Dictionary of pairs of hexagons with their distances only containing the ones that are within the allowed distance
    """
    
    def find_h3_line(hex_start, hex_end, max_iterations=10):
        """Finds an H3 line between two hexagons, potentially using midpoints if necessary."""
        try:
            # Try to create a direct line first
            return h3.grid_path_cells(hex_start, hex_end)
        except:
            if max_iterations <= 0:
                return None  # Failed to find a line after maximum iterations
            
            # Calculate the midpoint between the centers of the two hexagons
            center_start = h3.cell_to_latlng(hex_start)
            center_end = h3.cell_to_latlng(hex_end)
            midpoint = [(center_start[0] + center_end[0]) / 2, (center_start[1] + center_end[1]) / 2]
            
            # Convert midpoint to the nearest H3 hexagon
            midpoint_hex = h3.latlng_to_cell(midpoint[0], midpoint[1], h3.get_resolution(hex_start))
            
            # Recursively attempt to find lines using the midpoint
            first_half = find_h3_line(hex_start, midpoint_hex, max_iterations - 1)
            second_half = find_h3_line(midpoint_hex, hex_end, max_iterations - 1)
            
            if first_half is not None and second_half is not None:
                # Combine the two halves, removing the duplicate midpoint hex
                return first_half[:-1] + second_half
            else:
                return None

    # Dictionaries to save barrier lines, hexagons with distances, and direct neighbor distances
    barrier_lines = defaultdict(float)
    barrier_hex = defaultdict(list)
    hex_dist_to_direct_neighbors = defaultdict(list)
    new_time_bin = defaultdict(float)
    
    # Loop over all pairs of hexagons in the time bin
    for pair, distance in time_bin.items():
        pair = list(pair)
        # Check if the pair are direct neighbors
        if pair[0] in h3.grid_ring(pair[1], 1)[1]:
            # Get the line between the two hexagons
            boundary1 = h3.cell_to_boundary(pair[0])
            boundary2 = h3.cell_to_boundary(pair[1])
            # Get the pair of dots that the two hexagons share
            shared_boundary = frozenset([x for x in boundary1 if x in boundary2])
            # Add the line and its distance to the dictionary
            barrier_lines[shared_boundary] = distance
            hex_dist_to_direct_neighbors[pair[0]].append(distance)
            hex_dist_to_direct_neighbors[pair[1]].append(distance)
            # add the pair and the distance to the new time bin
            new_time_bin[frozenset(pair)] = distance
        else:
            # Get the line between the two hexagons using the find_h3_line function
            line = find_h3_line(pair[0], pair[1])
            if line is not None and len(line) <= allowed_distance:
                for hex in line:
                    barrier_hex[hex].append(distance)
                # add the pair and the distance to the new time bin
                new_time_bin[frozenset(pair)] = distance

    # Calculate the average distance for each hexagon and round it to 2 decimal places
    barrier_hex = {hex: round(sum(distances) / len(distances), 2) for hex, distances in barrier_hex.items()}
    
    # Get all distances for every hexagon that are not yet in the hex_dist_to_direct_neighbors to check for isolated hexagons
    for pair, distance in time_bin.items():
        for hex in pair:
            hex_dist_to_direct_neighbors[hex].append(distance)
    
    # extract the hexagons that are isolated given the threshold
    isolated_hex = [hex for hex, distances in hex_dist_to_direct_neighbors.items() if all(d >= threshold for d in distances)]
    
    return isolated_hex, barrier_lines, barrier_hex, new_time_bin


def find_closest_population(df, time_bin_index, isolated_hex, dist_matrix, threshold, gen_distances_pred, resolution):
    """
    Find the closest population for each isolated hexagon.

    Parameters:
    - df (pd.DataFrame): DataFrame containing the data with hexagon IDs and age groups.
    - time_bin_index (int): Index of the time bin to be analyzed.
    - isolated_hex (list): List of isolated hexagons.
    - dist_matrix (pd.DataFrame): A DataFrame representing the distance matrix.
    - threshold (float): Distance threshold to consider a hexagon isolated.
    - gen_distance_pred (np.array): Array of predicted genetic distances based on geographic distances.

    Returns:
    - tuple:
        - closest_populations (dict): Dictionary where keys are tuples of (isolated_hex, closest_hex) and values are distances.
        - new_isolated_hex (list): List of isolated hexagons that have no close population.
    """

    # Convert 'AgeGroup' column values to tuples of integers representing the start and end years
    df['AgeGroupTuple'] = df['AgeGroup'].apply(lambda x: tuple(map(int, x.split('-'))))
    
    # Sort the unique age group tuples to process them in chronological order
    time_bins = sorted(df['AgeGroupTuple'].unique())
    
    # Get the samples in the time bin of interest
    time_bin_df = df[df['AgeGroupTuple'] == time_bins[time_bin_index]]
    
    # Get column name for hexagons (it should be the only column with 'hex' in the name)
    hex_col = time_bin_df.filter(like='hex').columns[0]
    
    # Get all unique hexagons from the DataFrame
    hexagons = time_bin_df[hex_col].unique()
    
    # Get the samples in each hexagon
    samples_in_hex = time_bin_df.groupby(hex_col)['ID'].apply(list).to_dict()
    
    # Create a submatrix of the distance matrix for the samples in the hexagons
    all_samples = [sample for samples in samples_in_hex.values() for sample in samples]
    dist_matrix = dist_matrix.loc[all_samples, all_samples]
    
    # Prepare dictionary to hold the distances between the hexagons
    closest_populations = {}
    
    # List to hold the isolated hexagons that have no close population
    new_isolated_hex = []
    
    # Loop over all isolated hexagons
    for iso in isolated_hex:
        # Initialize the minimum distance to the threshold
        min_dist = float('inf')
        all_dist = {}
        closest_hex = None
        
        # get distances to every hexagon in the time bin
        Ids_in_iso = samples_in_hex.get(iso, [])
        for hex in hexagons:
            if hex == iso:
                continue

            Ids_in_neighbor = samples_in_hex.get(hex, [])
            distance = calc_avg_dist(Ids_in_iso, Ids_in_neighbor, dist_matrix)
            
            if distance < min_dist:
                min_dist = distance
                closest_hex = hex
            
        # add the distance to the dictionary with the two hexagons as a pair
        all_dist[frozenset([iso, closest_hex])] = round(distance, 2)
        
        # scale the distances by the estimated genetic differences
        scaled_distances = scale_distances(all_dist, gen_distances_pred, resolution)[0]
        
        if scaled_distances[frozenset([iso, closest_hex])] < threshold:
            pair = (iso, closest_hex)
            closest_populations[pair] = round(min_dist, 2)
        else:
            new_isolated_hex.append(iso)
    
    return closest_populations, new_isolated_hex


def impute_missing_hexagons(barrier_hex, num_runs=5):
    """
    Impute missing hexagons by iteratively adding neighboring hexagons that meet criteria.

    Parameters:
    - barrier_hex (dict): Dictionary with hexagon IDs as keys and average distances as values.
    - num_runs (int): Number of iterations to perform the imputation.

    Returns:
    - imputed_hex (dict): Dictionary with newly imputed hexagons and their average distances.
    """
    def impute(barrier_hex):
        """
        Impute missing hexagons for a single run based on neighbors.

        Parameters:
        - barrier_hex (dict): Dictionary with hexagon IDs as keys and average distances as values.

        Returns:
        - output_hex (dict): Dictionary of newly imputed hexagons with calculated average distances.
        """
        # Convert barrier_hex keys to a set for faster lookup
        barrier_hex_set = set(barrier_hex)
        new_barrier_hex = defaultdict(list)

        # Iterate through all barrier hexagons
        for hexagon in barrier_hex:
            # Find neighbors that are not in the barrier_hex
            neighbors = [hex for hex in h3.grid_disk(hexagon, 1) if hex not in barrier_hex_set]
            # Collect distances for these neighbors
            for neighbor in neighbors:
                new_barrier_hex[neighbor].append(barrier_hex[hexagon])

        # Calculate average distance for neighbors with at least 3 entries
        output_hex = {hexagon: round(sum(distances) / len(distances), 2)
                        for hexagon, distances in new_barrier_hex.items() if len(distances) >= 3}
        return output_hex

    # Create a copy of the barrier_hex
    imputed_hex = barrier_hex.copy()

    # Perform imputation for the specified number of runs
    for _ in range(num_runs):
        imputed_hex.update(impute(imputed_hex))

    # Remove hexagons that were originally in the barrier_hex
    imputed_hex = {hex: dist for hex, dist in imputed_hex.items() if hex not in barrier_hex}

    return imputed_hex


def scale_distances(time_bin, exsiting_pred=None, resolution=3):
    """
    Scales genetic distances by the estimated genetic differences modeled by a LOESS function of the geographic distances.

    Parameters:
    - time_bin: A dictionary where keys are pairs of hexagons and values are genetic distances between them.

    Returns:
    - output: A dictionary where keys are pairs of hexagons and values are scaled genetic distances.
    """

    def get_km_distance(hex1, hex2=None, resolution=3):
        """
        Calculates the geographic distance in kilometers between two hexagons.

        Parameters:
        - hex1: H3 index of the first hexagon.
        - hex2: H3 index of the second hexagon.

        Returns:
        - The distance in kilometers between the two hexagons.
        """
        # if it is the distance to itselve return a fiexed value based on the resoultion
        if hex2 is None:
            return 1281/(2.65**resolution)
        # else calculate the distance between the two hexagons based on the haversine formula
        coord1 = h3.cell_to_latlng(hex1)
        coord2 = h3.cell_to_latlng(hex2)
        return haversine(coord1, coord2)

    # Calculate km distances between the hexagons
    km_time_bin = {pair: get_km_distance(*pair) for pair in time_bin}

    # Convert genetic and geographic distances to numpy arrays
    gen_distances = np.array(list(time_bin.values()))
    geo_distances = np.array(list(km_time_bin.values()))
    # if there is no existing prediction, create one
    if exsiting_pred is None:
        # Apply LOESS smoothing to the genetic distances based on geographic distances
        lowess = sm.nonparametric.lowess
        gen_distances_pred = lowess(gen_distances, geo_distances, frac=0.5)
    else:
        gen_distances_pred = exsiting_pred

    # Scale genetic distances by the predicted values from the LOESS model
    output = {}
    for pair in time_bin:
        km_distance = km_time_bin[pair]
        # if km_distance not in gen_distances_pred[:, 0]: find the one closest to it
        if km_distance not in gen_distances_pred[:, 0]:
            km_distance = gen_distances_pred[:, 0][np.argmin(np.abs(gen_distances_pred[:, 0] - km_distance))]
        gen_distance = time_bin[pair]
        gen_distance_pred = gen_distances_pred[gen_distances_pred[:, 0] == km_distance][:, 1][0]
        # check if the predicted distance is 0
        if gen_distance == 0:
            output[pair] = 0
        # else get the log 2 of the ratio of the genetic distance to the predicted genetic distance
        else:
            output[pair] = round(math.log2(gen_distance / gen_distance_pred),2)

    return output, gen_distances_pred

def get_distance_lines(time_bin):
    """
    Converts all distances between hexagons to drawable lines.

    Parameters:
    - time_bin: A dictionary where keys are pairs of hexagons and values are distances between them.

    Returns:
    - lines: A dictionary where keys are pairs of coordinates (as frozensets) representing lines between hexagons and values are the distances rounded to two decimal places.
    """
    lines = {}
    for key, value in time_bin.items():
        hex1, hex2 = key
        coord1 = h3.cell_to_latlng(hex1)
        coord2 = h3.cell_to_latlng(hex2)
        # Create a frozenset of coordinates to represent the line
        shared_boundary = frozenset([coord1, coord2])
        lines[shared_boundary] = value
    return lines


def rename_time_bins(df):
    """
    Renames the time bins in a DataFrame to a more readable format.

    Parameters:
    - df (pd.DataFrame): DataFrame containing an 'AgeGroup' column with time bins as strings.

    Returns:
    - renamed_bins (list): List of renamed time bins in a more readable format.
    """
    # Convert the 'AgeGroup' column values to tuples of integers representing the start and end years
    df['AgeGroupTuple'] = df['AgeGroup'].apply(lambda x: tuple(map(int, x.split('-'))))
    
    # Sort the unique age group tuples to process them in chronological order
    time_bins = sorted(df['AgeGroupTuple'].unique())
    
    # Rename the time bins using the rename_times function
    renamed_bins = [rename_times(time_bin) for time_bin in time_bins]
    
    return renamed_bins


def rename_times(time_bin):
    """
    Renames a time bin into a more readable format, converting years to BC/AD notation.

    Parameters:
    - time_bin (tuple): A tuple of integers representing the start and end years.

    Returns:
    - str: A string representing the time bin in a more readable BC/AD format.
    """
    renamed_years = []
    for year in time_bin:
        if year < 1950:
            year = 1950 - year
            year_str = f"{year} AD"
        else:
            year = year - 1950
            year_str = f"{year} BC"
        renamed_years.append(year_str)
    
    return " - ".join(renamed_years)

def create_equal_age_groups(df, number_of_bins):
    """
    Creates n time bins with equally sized age groups.
    
    Parameters:
    - df (pd.DataFrame): DataFrame containing sample data with an 'Age' column.
    - number_of_bins (int): Number of bins to create.
    
    Returns:
    - age_groups (list): List of DataFrames, each representing a time bin.
    """
    temp_df = df.sort_values('Age')
    # Leave out the oldest 144 samples for creation of time bin size so we have a time range of 14000 years
    temp_df = temp_df.iloc[:-144]
    min_age = temp_df['Age'].min()
    max_age = temp_df['Age'].max()
    age_range = max_age - min_age
    bin_size = int(age_range / number_of_bins)
    age_groups = []
    low_b = min_age
    up_b = min_age + bin_size

    for _ in range(number_of_bins - 1):
        temp_df = df[(df['Age'] >= low_b) & (df['Age'] < up_b)]
        while temp_df.shape[0] < 5:
            up_b += 500
            temp_df = df[(df['Age'] >= low_b) & (df['Age'] < up_b)]
        age_groups.append(temp_df)
        low_b = up_b
        up_b += bin_size

    temp_df = df[(df['Age'] >= low_b)]
    age_groups.append(temp_df)

    return age_groups


def create_age_groups(df, number_of_bins):
    """
    Creates n time bins with equally distributed number of samples in each group.
    
    Parameters:
    - df (pd.DataFrame): DataFrame containing sample data with an 'Age' column.
    - number_of_bins (int): Number of bins to create.
    
    Returns:
    - age_groups (list): List of DataFrames, each representing a time bin.
    """
    total_samples = df.shape[0]
    sample_per_bin, remainder = divmod(total_samples, number_of_bins)
    temp_df = df.sort_values('Age')
    age_groups = []
    start_index = 0

    for bin_num in range(number_of_bins):
        end_index = start_index + sample_per_bin + (1 if remainder > 0 else 0)
        if remainder > 0:
            remainder -= 1
        age_groups.append(temp_df.iloc[start_index:end_index])
        start_index = end_index

    return age_groups


def name_age_groups(age_groups):
    """
    Creates a dictionary mapping sample IDs to age group names.
    
    Parameters:
    - age_groups (list): List of DataFrames, each representing a time bin.
    
    Returns:
    - name_dict (dict): Dictionary mapping sample IDs to age group names.
    """
    name_dict = {}
    for group in age_groups:
        min_age = group['Age'].min()
        max_age = group['Age'].max()
        for i in range(len(group)):
            id = group.iloc[i]['ID']
            name_dict[id] = f"{min_age}-{max_age}"
    return name_dict


def add_age_group_column(df, name_dict):
    """
    Adds an age group column to the DataFrame based on the name dictionary.
    
    Parameters:
    - df (pd.DataFrame): DataFrame containing sample data.
    - name_dict (dict): Dictionary mapping sample IDs to age group names.
    
    Returns:
    - df (pd.DataFrame): DataFrame with an added 'AgeGroup' column.
    """
    df['AgeGroup'] = df['ID'].map(name_dict)
    return df


def filter_df(df):
    """
    Filters the DataFrame for errors in the latitude and longitude columns.
    
    Parameters:
    - df (pd.DataFrame): DataFrame containing sample data with 'Latitude' and 'Longitude' columns.
    
    Returns:
    - filtered_df (pd.DataFrame): Filtered DataFrame.
    """
    filtered_df = df[df['Latitude'] != '..']
    filtered_df = filtered_df[filtered_df['Longitude'] != '..']
    return filtered_df


def assign_hexagon_to_samples(df, resolution):
    """
    Assigns a hexagon to each sample based on latitude and longitude.
    
    Parameters:
    - df (pd.DataFrame): DataFrame containing sample data with 'Latitude' and 'Longitude' columns.
    - resolution (int): H3 resolution for hexagon assignment.
    
    Returns:
    - df (pd.DataFrame): DataFrame with an added hexagon column.
    """
    hex_col = 'hex_res_' + str(resolution)
    df[hex_col] = df.apply(lambda x: h3.latlng_to_cell(float(x['Latitude']), float(x['Longitude']), res=resolution), axis=1)
    return df


def write_df(df, path):
    """
    Writes the DataFrame to a file.
    
    Parameters:
    - df (pd.DataFrame): DataFrame to be written.
    - path (str): File path for output.
    
    Returns:
    - None
    """
    df.to_csv(path, sep="\t", index=False)
    
    
def get_samples_per_time_bin(df):
    """
    Create a new dataframe with the number of samples per time bin.
    
    Parameters:
    - df (pd.DataFrame): DataFrame containing the data.
    
    Returns:
    - new_df (pd.DataFrame): DataFrame with the number of samples per time bin.
    """
    # Create a new DataFrame with the number of samples per time bin
    new_df = df.groupby('AgeGroupTuple')['ID'].count().reset_index()
    
    # Rename the columns
    new_df.columns = ['Time Bin', 'Number of samples']
    
    # bring the AgeGroup column in to a more readable format
    new_df['Time Bin'] = new_df['Time Bin'].apply(rename_times)
    
    return new_df


def label_samples(path, number_of_bins=20, resolution=2, equally_sized=False):
    """
    Labels the ancient samples with the time bins and hexagons and returns the DataFrame.
    
    Parameters:
    - path (str): Path to the input file.
    - number_of_bins (int): Number of time bins to create. Default is 20.
    - resolution (int): H3 resolution for hexagon assignment. Default is 2.
    - equally_sized (bool): Whether to create equally sized age groups. Default is False.
    
    Returns:
    - new_df (pd.DataFrame): DataFrame with added age group and hexagon columns.
    """
    df = read_df(f'{path}/0_data/Ancient_samples.txt')
    
    # Create the age groups
    if equally_sized:
        age_groups = create_equal_age_groups(df, number_of_bins)
    else:
        age_groups = create_age_groups(df, number_of_bins)
    
    # Get the name for each age group
    name_dict = name_age_groups(age_groups)
    
    # Create a new DataFrame with the age group column
    new_df = add_age_group_column(df, name_dict)
    
    # Filter for errors in the latitude and longitude columns
    new_df = filter_df(new_df)
    
    # Assign a hexagon to each sample
    new_df = assign_hexagon_to_samples(new_df, resolution=resolution)
    
    # Write the DataFrame to a file
    write_df(new_df, f'{path}/0_data/Ancient_samples_with_time_hexagon.txt')
    
    return new_df