# This script reads the ancient samples with labels and divides them into age groups.
# The age groups are created by dividing the samples into n bins with an equal number of samples in each group.
# The age group is then added to the data frame and written to a file.

import pandas as pd
from h3 import h3
from func import read_df

# function that creates n time bins with equally sized age groups
def create_equal_age_groups(df, number_of_bins):
    temp_df = df.sort_values('Age')
    # leave out the oldest 60 samples for creation of time bin size
    temp_df = temp_df.iloc[:-60]
    # get the smallest and largest age
    min_age = temp_df['Age'].min()
    max_age = temp_df['Age'].max()

    # calculate the range of ages
    age_range = max_age - min_age
    # calculate the size of each bin
    bin_size = int(age_range/number_of_bins)
    age_groups = [] 
    low_b = min_age
    up_b = min_age + bin_size

    for _ in range(number_of_bins-1):
        # Use a single step to filter the DataFrame for the age range
        temp_df = df[(df['Age'] >= low_b) & (df['Age'] < up_b)]
        age_groups.append(temp_df)  # Append the filtered DataFrame to your list
        
        # Update bounds for the next iteration
        low_b += bin_size
        up_b += bin_size
        
    temp_df = df[(df['Age'] >= low_b)]
    age_groups.append(temp_df)
        
    return age_groups


# function that creates n time bins with equally distributed number of samples in each group
def create_age_groups(df, number_of_bins):
    total_samples = df.shape[0]
    sample_per_bin, remainder = divmod(total_samples, number_of_bins)
    temp_df = df.sort_values('Age')
    age_groups = []
    start_index = 0
    # iterate over the number of bins
    for bin_num in range(number_of_bins):
        # Calculate the end index for the current bin; add an extra sample to some bins to distribute the remainder
        if remainder > 0:
            end_index = start_index + sample_per_bin + 1
            remainder -= 1
        else:
            end_index = start_index + sample_per_bin
        
        # Append the subset of the dataframe corresponding to the current bin
        age_groups.append(temp_df.iloc[start_index:end_index])
        
        # Update the start index for the next bin
        start_index = end_index
        
    return age_groups

# function that creates a name for each age group
def name_age_groups(age_groups):
    name_dict = {}
    for group in age_groups:
        min_age = group['Age'].min()
        max_age = group['Age'].max()
        for i in range(len(group)):
            id = group.iloc[i]['ID']
            name_dict[id] = f"{min_age}-{max_age}"
    return name_dict

# function that adds the age group to the dataframe
def add_age_group_column(df, name_dict):
    df['AgeGroup'] = df['ID'].map(name_dict)
    return df

# function that filters the data frame for errors in the latitude and longitude columns
def filter_df(df):
    # filter out the samples with missing latitude and longitude
    filtered_df = df[df['Latitude'] != '..']
    filtered_df = filtered_df[filtered_df['Longitude'] != '..']
    return filtered_df

# function that assigns a hexagon to each sample with the given resolution
def assign_hexagon_to_samples(df, resolution):
    hex_col = 'hex_res_'+str(resolution)
    df[hex_col] = df.apply(lambda x: h3.geo_to_h3(float(x['Latitude']), float(x['Longitude']), resolution=resolution), axis=1)
    return df

# function that writes the age groups to a file
def write_df(df, path):
    df.to_csv(path, sep="\t", index=False)
    return

# function that labels the ancient samples with the time bins and hexagons and returns the data frame
def label_samples(path, number_of_bins=30,resolution=2, equally_sized = False):
    df = read_df(f'{path}/0_data/Ancient_samples.txt')
    # create the age groups
    if equally_sized:
        age_groups = create_equal_age_groups(df, number_of_bins)
    else:
        age_groups = create_age_groups(df, number_of_bins)
    # get the name for each age group
    name_dict = name_age_groups(age_groups)
    # create a new dataframe with the age group column
    new_df = add_age_group_column(df, name_dict)
    # filter for errors in the latitude and longitude coulmns
    new_df= filter_df(new_df)
    # assign a hexagon to each sample
    new_df = assign_hexagon_to_samples(new_df, resolution=resolution)
    # write the dataframe to a file
    write_df(new_df, f'{path}/0_data/Ancient_samples_with_time_hexagon.txt')
    return new_df