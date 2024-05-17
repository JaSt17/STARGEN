# This script reads the ancient samples with labels and divides them into age groups.
# The age groups are created by dividing the samples into n bins with an equal number of samples in each group.
# The age group is then added to the data frame and written to a file.

import pandas as pd
from h3 import h3
from func import read_df


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
    # Leave out the oldest 100 samples for creation of time bin size
    temp_df = temp_df.iloc[:-100]
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
            up_b += 50
            temp_df = df[(df['Age'] >= low_b) & (df['Age'] < up_b)]
        age_groups.append(temp_df)
        low_b += bin_size
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
    df[hex_col] = df.apply(lambda x: h3.geo_to_h3(float(x['Latitude']), float(x['Longitude']), resolution=resolution), axis=1)
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