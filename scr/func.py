import pandas as pd
import numpy as np
from h3 import h3
from collections import defaultdict

# function that reads the ancient DNA Annotations file into a data frame
def read_df(path):
    df = pd.read_csv(path, sep="\t")
    #convert the age column to numeric
    df['Age'] = pd.to_numeric(df['Age'], errors='coerce')
    return df


# function that calculates the average distance between two groups of samples
def calc_avg_dist(samples_hex1, samples_hex2, dist_matrix):
    return dist_matrix.loc[samples_hex1, samples_hex2].values.flatten().mean()


def calc_neighbor_dist(hexagons, dist_matrix, time_bin_df, hex_col, k_neighbors = 1, scale_by_distance=False):
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
        if [len(neighbors[k]) for k in neighbors].count(0) == len(neighbors):
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
def calc_dist_time_bin(df, dist_matrix=None, k_neighbors=1, scale_by_distance=False):
    
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
        average_distances = calc_neighbor_dist(hexagons, dist_matrix, time_bin_df, hex_col, k_neighbors, scale_by_distance)

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


# function that gets isolated hexagons and barriers for each time bin
def get_isolated_hex_and_barriers(time_bin, hexagons, n, threshold):
    # directory to save the barrier lines and hexagons with their distances
    barrier_hex = defaultdict(list)
    barrier_lines = {}
    # dictionary to save hexagons and their direct neighbor distances to check if their was a migration barrier
    hex_dist_to_direct_neighbors = defaultdict(list)
    # loop over all pairs of hexagons in the time bin
    for pair in time_bin:
        distance = time_bin[pair]
        pair = list(pair)
        # check if the pair are direct neighbors
        if pair[0] in h3.k_ring_distances(pair[1], 1)[1]:
            # get the line between the two hexagons
            boundary1 = h3.h3_to_geo_boundary(pair[0])
            boundary2 = h3.h3_to_geo_boundary(pair[1])
            # get the pair of dots that the two hexagons share
            shared_boundary = frozenset([x for x in boundary1 if x in boundary2])
            # add the line and its distance to the dictionary
            barrier_lines[shared_boundary] = round(distance,2)
            hex_dist_to_direct_neighbors[pair[0]].append(distance)
            hex_dist_to_direct_neighbors[pair[1]].append(distance)
        # if the hexagons are further appart
        else:
            # try to draw a line between the two hexagons
            try:
                # if found add the hexagons to the barrier_hex dictionary and add the distance to the list
                line = h3.h3_line(pair[0], pair[1])[1:-1]
                for hex in line:
                    if hex not in hexagons:
                        # if the hexagon is already in the dictionary add the distance with a frequency of n-len(line) which accounts for the distance between the hexagons
                        dist_array = [distance] * max(n-len(line),1)
                        if hex in hexagons:
                            barrier_hex[hex].extend(dist_array)
                        else:
                            barrier_hex[hex] = dist_array
            except:
                continue

    # Calculate the average distance for each hexagon and round it to 2 decimal places
    barrier_hex = {hex: round(sum(distances)/len(distances), 2) for hex, distances in barrier_hex.items()}
        
    # Create a list of isolated hexagons
    # add isolated hexagons that have direct neighbors
    isolated_hex = [hex for hex, distances in hex_dist_to_direct_neighbors.items() if all(x >= threshold for x in distances)]
    # add isolated hexagons that have no direct neighbors
    isolated_hex += [hex for hex in hexagons if hex not in hex_dist_to_direct_neighbors and all(barrier_hex[n] >= threshold for n in h3.k_ring_distances(hex, 1)[1] if n in barrier_hex)]
    
    return isolated_hex, barrier_lines, barrier_hex


# function that finds the closest population for each isolated hexagon
def find_closest_population(df, time_bin, isolated_hex, dist_matrix, threshold):
    # Convert the 'AgeGroup' column values to tuples of integers representing the start and end years,
    df['AgeGroupTuple'] = df['AgeGroup'].apply(lambda x: tuple(map(int, x.split('-'))))
    # Sort the unique age group tuples to process them in a chronological order.
    time_bins = sorted(df['AgeGroupTuple'].unique())
    # get the samples in the time bin of interest
    time_bin_df = df[df['AgeGroupTuple'] == time_bins[time_bin]]
    # get column name for the hexagons (it should be the only column with 'hex' in the name)
    hex_col = time_bin_df.filter(like='hex').columns[0]
    # get all unique hexagons from the dataframe
    hexagons = time_bin_df[hex_col].unique()
    # get the samples in each hexagon
    samples_in_hex = time_bin_df.groupby(hex_col)['ID'].apply(list).to_dict()
    all_samples = [sample for samples in samples_in_hex.values() for sample in samples]
    # create a submatrix of the distance matrix for the samples in the hexagons
    dist_matrix = dist_matrix.loc[all_samples, all_samples]
    # empty dictrionary to hold the distances between the hexagons
    closest_populations = {}
    # empty list to hold the isolated hexagons that have no migration
    new_isolated_hex = []
    
    # loop over all isolated hexagons
    for iso in isolated_hex:
        # reset the min_dist and closest_hex
        closest_hex = None
        min_dist = threshold
        # check the distance to every hexagon in that time bin
        for hex in hexagons:
            # skip if the hexagon is the same as the isolated hexagon
            if hex == iso:
                continue
            
            Ids_in_hexagon = samples_in_hex.get(iso, [])
            Ids_in_neighbor = samples_in_hex.get(hex, [])
            # calculate the average distance between the hexagon and its neighbor
            distance = calc_avg_dist(Ids_in_hexagon, Ids_in_neighbor, dist_matrix)
            # check if the distance is lower than the current minimum distance
            if distance < min_dist:
                min_dist = distance
                closest_hex = hex
        # Add closest hexagon
        if closest_hex is not None:
            pair = frozenset([iso, closest_hex])
            closest_populations[pair] = round(min_dist, 2)
        else:
            new_isolated_hex.append(iso)
    
    return closest_populations, new_isolated_hex


def impute_missing_hexagons_multiple_runs(barrier_hex, hexagons, num_runs=5):
    # impute the missing hexagons in the barrier_hex
    def impute_missing_hexagons(barrier_hex):
        # convert the barrier_hex to sets for faster lookup
        barrier_hex_set = set(barrier_hex)
        # create a dictionary to store
        new_barrier_hex = defaultdict(list)
        # loop through all the barrier hexagons
        for hexagon in barrier_hex:
            # check if the neighbors are not in the barrier_hex
            neighbors = [hex for hex in h3.k_ring(hexagon, 1) if hex not in barrier_hex_set]
            # loop over all remaining neighbors
            for neighbor in neighbors:
                if neighbor in new_barrier_hex:
                    new_barrier_hex[neighbor].append(barrier_hex[hexagon])
                else:
                    new_barrier_hex[neighbor] = [barrier_hex[hexagon]]
        output_hex = defaultdict(float)
        # delete the hexagons that have less than 3 neighbors
        for hexagon, distances in new_barrier_hex.items():
            if len(distances) < 3:
                continue
            else:
                output_hex[hexagon] = round(sum(distances)/len(distances), 2)
        return output_hex

    # create a copy of the barrier_hex
    imputed_hex = barrier_hex.copy()
    # loop through the number of runs
    for _ in range(num_runs):
        imputed_hex.update(impute_missing_hexagons(imputed_hex))
        
    # delete the hexagons that are already in the barrier_hex
    imputed_hex = {hex: dist for hex, dist in imputed_hex.items() if hex not in barrier_hex}

    return imputed_hex

        
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