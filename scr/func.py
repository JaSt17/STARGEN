import pandas as pd
import numpy as np
from h3 import h3

# function that reads the ancient DNA Annotations file into a data frame
def read_df(path):
    df = pd.read_csv(path, sep="\t")
    #convert the age column to numeric
    df['Age'] = pd.to_numeric(df['Age'], errors='coerce')
    return df

# function that calculates the average distance between two groups of samples
def calc_avg_dist(samples_hex1, samples_hex2, dist_matrix):
    return np.mean(dist_matrix.loc[samples_hex1, samples_hex2].values.flatten())

def calc_neighbor_dist(hexagons, dist_matrix, time_bin_df, hex_col, k_neighbors = 1, allow_k_distance=False, scale_by_distance=False):
    # get the samples in each hexagon
    samples_in_hex = time_bin_df.groupby(hex_col)['ID'].apply(list).to_dict()
    # create a list of all values in samples_in_hex
    all_samples = [sample for samples in samples_in_hex.values() for sample in samples]
    # create a submatrix of the distance matrix for the samples in the hexagons
    dist_matrix = dist_matrix.loc[all_samples, all_samples]
    # initialize the dictionary to store the average distances between neighboring hexagons
    averages = {}
    # get the set of hexagons
    hexagons_set = set(hexagons)
    # initialize the cache for the average distances
    avg_dist_cache = {}
    # initialize the cache for the k-ring distances
    k_ring_distances_cache = {}

    for hexagon in hexagons:
        neighbors = dict()
        # get the neighbors of the hexagon in k distance
        for k in range(1, k_neighbors+1):
            # if the neighbors have not been calculated yet, calculate them
            if (hexagon, k) not in k_ring_distances_cache:
                k_ring_distances_cache[(hexagon, k)] = set(h for h in h3.k_ring_distances(hexagon, k)[k] if h in hexagons_set)
            neighbors[k] = k_ring_distances_cache[(hexagon, k)]
        # if there are no neighbors in k distance, and the user allows for more than k distance, get the neighbors in 20 distance
        if allow_k_distance and [len(neighbors[k]) for k in neighbors].count(0) == len(neighbors):
            k = k_neighbors + 1
            while all(len(neighbors[k]) == 0 for k in neighbors) and k < 20:
                # if the neighbors have not been calculated yet, calculate them
                if (hexagon, k) not in k_ring_distances_cache:
                    k_ring_distances_cache[(hexagon, k)] = set(h for h in h3.k_ring_distances(hexagon, k)[k] if h in hexagons_set)
                neighbors[k] = k_ring_distances_cache[(hexagon, k)]
                k += 1
                
        # calculate the average distance between the hexagon and its neighbors
        for k in neighbors.keys():
            for neighbor in neighbors[k]:
                Ids_in_hexagon = samples_in_hex.get(hexagon, [])
                Ids_in_neighbor = samples_in_hex.get(neighbor, [])
                # get the pair of hexagons
                pair = frozenset([hexagon, neighbor])

                # check if the average distance has already been calculated
                if pair not in avg_dist_cache:
                    # calculate the average distance between the hexagon and its neighbor
                    distance = calc_avg_dist(Ids_in_hexagon, Ids_in_neighbor, dist_matrix)
                    # scale the distance by the distance between the hexagon and its neighbor
                    if scale_by_distance:
                        distance = distance / (0.9 + k/10)
                    avg_dist_cache[pair] = distance
                    
                averages[pair] = round(avg_dist_cache[pair], 2)

    return averages

# this function calculates the average distance between the each hexagon and its neighbors for each time bin
def calc_dist_time_bin(df, dist_matrix=None, k_neighbors=1, allow_k_distance=False, scale_by_distance=False):
    
    # get column name for the hexagons (it should be the only column with 'hex' in the name)
    hex_col = str(df.columns[df.columns.str.contains('hex')][0])
    
    # Convert the 'AgeGroup' column values to tuples of integers representing the start and end years,
    df['AgeGroupTuple'] = df['AgeGroup'].apply(lambda x: tuple(map(int, x.split('-'))))
    
    # Sort the unique age group tuples to process them in a chronological order.
    time_bins = sorted(df['AgeGroupTuple'].unique())
    averages = {}
    # Iterate over each time bin.
    for time_bin in time_bins:
        # Format the current time bin as a string for labeling purposes.
        bin_label = rename_times(time_bin)
        
        # get subset of the data frame for that time bin
        time_bin_df = df[df['AgeGroupTuple'] == time_bin]

        # get all unique hexagons for that time bin
        hexagons = time_bin_df[hex_col].unique()
        
        # Calculate the average distance for each hexagon to its neighbors within the current time bin.
        average_distances = calc_neighbor_dist(hexagons, dist_matrix, time_bin_df, hex_col, k_neighbors, allow_k_distance, scale_by_distance)

        # Append the calculated average distances to the dictionary, using the time bin label as the key.
        averages.update({bin_label: average_distances})

    # Return the dictionary with the average distances between neighboring hexagons for each time bin.
    return averages

# function that gets the hexagons for each time bin
def get_time_bin_hexagons(df):
    
    # get column name for the hexagons (it should be the only column with 'hex' in the name)
    hex_col = str(df.columns[df.columns.str.contains('hex')][0])
    
    # Convert the 'AgeGroup' column values to tuples of integers representing the start and end years,
    df['AgeGroupTuple'] = df['AgeGroup'].apply(lambda x: tuple(map(int, x.split('-'))))
    
    # Sort the unique age group tuples to process them in a chronological order.
    time_bins = sorted(df['AgeGroupTuple'].unique())
    
    time_bin_hexagons = {}
    # Iterate over each time bin.
    for time_bin in time_bins:
        # Format the current time bin as a string for labeling purposes.
        bin_label = rename_times(time_bin)
        
        # get data frame for that time bin
        time_bin_df = df[df['AgeGroupTuple'] == time_bin]

        # get all unique hexagons for that time bin
        hexagons = time_bin_df[hex_col].unique().tolist()
        
        # add the hexagons to the dictionary
        time_bin_hexagons.update({bin_label: hexagons})
        
    return time_bin_hexagons
        
# function that normalizes the distance values on a interval from 0 to 1
def normalize_distances(time_bin):
    normalized_time_bin = time_bin.copy()
    min_dist, max_dist = get_min_max_dist(time_bin)
    for pair in time_bin:
        normalized_time_bin[pair] = round((time_bin[pair] - min_dist) / (max_dist - min_dist), 5)
    return normalized_time_bin

# function that gets the minimum and maximum distance values for a time bin
def get_min_max_dist(time_bin):
    min_dist = 1
    max_dist = 0
    for pair in time_bin:
        if time_bin[pair] < min_dist:
            min_dist = time_bin[pair]
        if time_bin[pair] > max_dist:
            max_dist = time_bin[pair]
    return min_dist, max_dist

# function that renames the time bins into a more readable format
def rename_time_bins(df):
    # Convert the 'AgeGroup' column values to tuples of integers representing the start and end years,
    df['AgeGroupTuple'] = df['AgeGroup'].apply(lambda x: tuple(map(int, x.split('-'))))
    # Sort the unique age group tuples to process them in a chronological order.
    time_bins = sorted(df['AgeGroupTuple'].unique())
    renamed_bins = []
    for time_bin in time_bins:
        renamed_bins.append(rename_times(time_bin))
    return renamed_bins

# function that renames a time bin into a more readable format
def rename_times(time_bin):
    renamed_years = []
    # get years from time_bin
    for year in time_bin:
        # the time in the dataset is measured from 1950
        if int(year) < 1950:
            year = 1950 - int(year)
            year = str(year) + " AD"
        else:
            year = int(year) - 1950
            year = str(year) + " BC"
        renamed_years.append(year)
    return(" - ".join(renamed_years))  # Append to renamed_bins