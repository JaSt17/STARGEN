import pandas as pd
import numpy as np
from h3 import h3
from collections import defaultdict
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
from pykrige.ok import OrdinaryKriging
from haversine import haversine
from scipy.optimize import curve_fit

# function that reads the ancient DNA Annotations file into a data frame
def read_df(path):
    df = pd.read_csv(path, sep="\t")
    #convert the age column to numeric
    df['Age'] = pd.to_numeric(df['Age'], errors='coerce')
    return df


# function that calculates the average distance between two groups of samples
def calc_avg_dist(samples_hex1, samples_hex2, dist_matrix):
    return dist_matrix.loc[samples_hex1, samples_hex2].values.flatten().mean()


def calc_neighbor_dist(hexagons, dist_matrix, time_bin_df, hex_col):
    # get the samples in each hexagon
    samples_in_hex = time_bin_df.groupby(hex_col)['ID'].apply(list).to_dict()
    # create a list of all values in samples_in_hex
    all_samples = [sample for samples in samples_in_hex.values() for sample in samples]
    # create a submatrix of the distance matrix for the samples in the hexagons
    dist_matrix = dist_matrix.loc[all_samples, all_samples]
    # initialize the dictionary to store the average distances between neighboring hexagons
    averages = {}
    
    # function which gets the neighbors for all hexagons usin Delauny triangulation
    def get_neighbors(hexagons):
        # get the centroid of every hexagon
        coord = []
        for hex in hexagons:
            coord.append(h3.h3_to_geo(hex))
        # change coord to and hexagons to np array
        coord = np.array(coord)
        hexagons = np.array(hexagons)
        # calc the Delaunay triangle
        tri = Delaunay(coord)
        # function to get the edges
        def tris2edges(tris):
            edges = set([])
            for tri in tris:
                for k in range(3):
                    i,j = tri[k], tri[(k+1)%3]
                    i,j = min(i,j), max(i,j)  # canonical edge representation
                    edges.add((i,j))
            return edges
        # get the edges (pairs of hexagons)
        edges = tris2edges(tri.simplices)
        # use the edges to create a dictionary with all the neighborhoods that we want to calculate
        neighbors = {}
        for (i,j) in edges:
            if hexagons[i] in neighbors:
                neighbors[hexagons[i]].append(hexagons[j])
            else:
                neighbors[hexagons[i]] = [hexagons[j]]
        return neighbors
    
    # get the neighbors for whom we will calculate the distances
    neighbors = get_neighbors(hexagons)
                
    # calculate the average distance between the hexagon and its neighbors
    for hexagon in neighbors.keys():
        for neighbor in neighbors[hexagon]:
            Ids_in_hexagon = samples_in_hex.get(hexagon, [])
            Ids_in_neighbor = samples_in_hex.get(neighbor, [])
            # get the pair of hexagons
            pair = frozenset([hexagon, neighbor])

            # calculate the average distance between the hexagon and its neighbor
            distance = calc_avg_dist(Ids_in_hexagon, Ids_in_neighbor, dist_matrix)
                
            averages[pair] = round(distance, 2)

    return averages


# this function calculates the average distance between the each hexagon and its neighbors for each time bin
def calc_dist_time_bin(df, dist_matrix=None):
    
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
        average_distances = calc_neighbor_dist(hexagons, dist_matrix, time_bin_df, hex_col)

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
def get_isolated_hex_and_barriers(time_bin, hexagons, threshold, allowed_distance=12):
    
    def find_h3_line(hex_start, hex_end, max_iterations=10):
        """Finds an H3 line between two hexagons, potentially using midpoints if necessary."""
        try:
            # Try to create a direct line first
            return h3.h3_line(hex_start, hex_end)
        except:
            if max_iterations <= 0:
                return None  # Failed to find a line after maximum iterations
            
            # Calculate the midpoint between the centers of the two hexagons
            center_start = h3.h3_to_geo(hex_start)
            center_end = h3.h3_to_geo(hex_end)
            midpoint = [(center_start[0] + center_end[0]) / 2, (center_start[1] + center_end[1]) / 2]
            
            # Convert midpoint to the nearest H3 hexagon
            midpoint_hex = h3.geo_to_h3(midpoint[0], midpoint[1], h3.h3_get_resolution(hex_start))
            
            # Recursively attempt to find lines using the midpoint
            first_half = find_h3_line(hex_start, midpoint_hex, max_iterations - 1)
            second_half = find_h3_line(midpoint_hex, hex_end, max_iterations - 1)
            
            if first_half is not None and second_half is not None:
                # Combine the two halves, removing the duplicate midpoint hex
                return first_half[:-1] + second_half
            else:
                return None

    # directory to save the barrier lines and hexagons with their distances
    barrier_hex = defaultdict(list)
    barrier_lines = defaultdict(float)
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
            # get the line between the two hexagons using the find_h3_line function
            line = find_h3_line(pair[0], pair[1])
            if line is not None:
                # check if the line is shorter than the allowed distance
                if len(line) <= allowed_distance:
                    for hex in line:
                        # if the hexagon is already in the dictionary, append the distance
                        if hex in barrier_hex:
                            barrier_hex[hex].append(distance)
                        # if not, add the hexagon and the distance
                        else:
                            barrier_hex[hex] = [distance]

    # Calculate the average distance for each hexagon and round it to 2 decimal places
    barrier_hex = {hex: round((sum(distances)/len(distances)), 2) for hex, distances in barrier_hex.items()}
    
    # get all distances for every hexagon
    for pair in time_bin:
        distance = time_bin[pair]
        for hex in pair:
            if hex not in hex_dist_to_direct_neighbors:
                hex_dist_to_direct_neighbors[hex] = [distance]
            else:
                hex_dist_to_direct_neighbors[hex].append(distance)
                
    isolated_hex = [hex for hex, distances in hex_dist_to_direct_neighbors.items() if all(x >= threshold for x in distances)]
    
    return isolated_hex, barrier_lines, barrier_hex


# function that finds the closest population for each isolated hexagon
def find_closest_population(df, time_bin, isolated_hex, dist_matrix, threshold):
    # Convert the 'AgeGroup' column values to tuples of integers representing the start and end years,
    df['AgeGroupTuple'] = df['AgeGroup'].apply(lambda x: tuple(map(int, x.split('-'))))
    # Sort the unique age group tuples to process them in a chronological order.
    time_bins = sorted(df['AgeGroupTuple'].unique())
    # Get the samples in the time bin of interest
    time_bin_df = df[df['AgeGroupTuple'] == time_bins[time_bin]]
    # Get column name for the hexagons (it should be the only column with 'hex' in the name)
    hex_col = time_bin_df.filter(like='hex').columns[0]
    # Get all unique hexagons from the dataframe
    hexagons = time_bin_df[hex_col].unique()
    # Get the samples in each hexagon
    samples_in_hex = time_bin_df.groupby(hex_col)['ID'].apply(list).to_dict()
    
    # Create a submatrix of the distance matrix for the samples in the hexagons
    all_samples = [item for sublist in samples_in_hex.values() for item in sublist]
    dist_matrix = dist_matrix.loc[all_samples, all_samples]
    
    # Prepare dictionary to hold the distances between the hexagons
    closest_populations = {}
    # List to hold the isolated hexagons that have no migration
    new_isolated_hex = []
    
    # Loop over all isolated hexagons
    for iso in isolated_hex:
        # Initialize the minimum distance to the threshold
        min_dist = threshold
        closest_hex = None
        
        # Calculate distances to every hexagon in the time bin
        Ids_in_iso = samples_in_hex.get(iso, [])
        for hex in hexagons:
            if hex == iso:
                continue

            Ids_in_neighbor = samples_in_hex.get(hex, [])
            distance = calc_avg_dist(Ids_in_iso, Ids_in_neighbor, dist_matrix)
            
            if distance < min_dist:
                min_dist = distance
                closest_hex = hex
        
        # Add the closest hexagon if found
        if closest_hex is not None:
            # Create a tuple with iso as the first element
            pair = (iso, closest_hex)
            # Add the pair to the dictionary with the distance
            closest_populations[pair] = round(min_dist, 2)
        else:
            new_isolated_hex.append(iso)
    
    return closest_populations, new_isolated_hex


# function that imputes missing hexagons using Kriging interpolation with a spherical variogram model
def impute_missing_hexagons(barrier_hex, num_runs=5):
    # Functions that are needed within the function:
    # Function to convert hexagon indices to lat/lon coordinates
    def hex_to_latlon(hex_index):
        return h3.h3_to_geo(hex_index)

    # Function to perform Kriging interpolation
    def kriging_interpolation(known_points, unknown_points):
        known_coords = np.array([hex_to_latlon(h) for h in known_points.keys()])
        known_values = np.array(list(known_points.values()))

        # Extract latitudes and longitudes separately for the Kriging function
        lats = known_coords[:, 0]
        lons = known_coords[:, 1]

        # Create OrdinaryKriging object
        OK = OrdinaryKriging(lons, lats, known_values, variogram_model='spherical', verbose=False, enable_plotting=False)

        # Interpolate unknown points
        unknown_coords = np.array([hex_to_latlon(h) for h in unknown_points])
        unknown_lats = unknown_coords[:, 0]
        unknown_lons = unknown_coords[:, 1]

        interpolated_values, _ = OK.execute('points', unknown_lons, unknown_lats)
        return dict(zip(unknown_points, np.round(interpolated_values, 2)))
    
    # Function to find unknown hexagons that are close to known hexagons
    def get_unknown_hexagons(barrier_hex, num_runs=10):
        known_hex_set = set(barrier_hex.keys())
        unknown_hex_set = set()

        for _ in range(num_runs):
            # Track all unknown hexagons to process
            unknown_hex_dict = {}
            
            for hexagon in known_hex_set:
                # Find neighbors not in known hex set
                neighbors = [hex for hex in h3.k_ring(hexagon, 1) if hex not in known_hex_set]
                for neighbor in neighbors:
                    if neighbor in unknown_hex_dict:
                        unknown_hex_dict[neighbor] += 1
                    else:
                        unknown_hex_dict[neighbor] = 1
            
            # keep only hexagons with at least 3 known neighbors
            unknown_hex_dict = {hex: count for hex, count in unknown_hex_dict.items() if count >= 3}
            # add the unknown hexagons to the sets
            known_hex_set.update(unknown_hex_dict.keys())
            unknown_hex_set.update(unknown_hex_dict.keys())
        
        return list(unknown_hex_set)
    
    # Main part of the function:
    new_hex = barrier_hex.copy()
    known_hex_set = set(new_hex.keys())
    unknown_hexagons = get_unknown_hexagons(barrier_hex, num_runs)
        
    imputed_hex = kriging_interpolation(new_hex, unknown_hexagons)
    
    return imputed_hex


# function scales genetic distances by the estimated genetic differents modeled by logistic growth of the geographic distances
def scale_distances(time_bin):
    # function that gets the km distance between two hexagons
    def get_km_distance(hex1, hex2):
            # get the centroid of the hexagons
            coord1 = h3.h3_to_geo(hex1)
            coord2 = h3.h3_to_geo(hex2)
            # get the distance between the two points
            return haversine(coord1, coord2)

    # get km distances between the hexagons
    km_time_bin = {pair: get_km_distance(*pair) for pair in time_bin}
    
    # normalize the km distances and genetic distances between 0 and 1 
    km_time_bin = normalize_distances(km_time_bin)
    time_bin = normalize_distances(time_bin)
    
    # get numpy arrays of genetic distances and km distances
    gen_distances = np.array(list(time_bin.values()))
    geo_distances = np.array(list(km_time_bin.values()))
    
    # Define the logistic growth function
    def logistic_growth(x, L, k, x_0):
        return L / (1 + np.exp(-k * (x - x_0)))

    # Initial guesses for the parameters
    initial_guesses = [max(gen_distances), 1, np.median(geo_distances)]

    # Fit the logistic growth function to the data
    params, covariance = curve_fit(logistic_growth, geo_distances, gen_distances, p0=initial_guesses)

    # Extract the fitted parameters
    L, k, x_0 = params
    
    # Use the predicted parameters to scale the genetic distances according to their geographic distance
    output = {}
    
    for pair in time_bin:
        hex1, hex2 = pair
        output[pair] = time_bin[pair]/logistic_growth(get_km_distance(hex1, hex2), *params)
    
    return output

# Function that converts all distances to drawable lines
def get_distance_lines(time_bin):
    lines = {}
    for key,value in time_bin.items():
        key = list(key)
        hex1 = key[0]
        hex2 = key[1]
        cord1 = h3.h3_to_geo(hex1)
        cord2 = h3.h3_to_geo(hex2)
        # get the pair of dots that the two hexagons share
        shared_boundary = frozenset([cord1, cord2])
        lines[shared_boundary] = round(value, 2)
    return lines

        
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