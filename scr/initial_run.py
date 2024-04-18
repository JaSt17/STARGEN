# This script is responsible for the initial run of the pipeline.
# It reads the excel file with the ancient samples, filters the data, and writes the ancient samples to a new file.
# It also calculates the euclidean distance matrix from the plink output and writes it to a pickle file.
# The script is called with the path to the excel file as an argument.

import pickle
import os
import sys
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform

# function that filters the data frame for errors in the latitude and longitude columns
def filter_df(df):
    # filter out the samples with missing latitude and longitude
    filtered_df = df[df['Lat.'] != '..']
    filtered_df = filtered_df[filtered_df['Long.'] != '..']
    # drop na values if there are any in Latitude and Longitude
    filtered_df = filtered_df.dropna(subset=['Lat.', 'Long.'])
    return filtered_df

# function that writes the ancient samples to a new file
def write_sample_list(df, path):
    # get the columns that are needed
    index = df['#']
    genetic_id = df['Genetic ID']
    countries_list = df['Political Entity']
    latitude = df['Lat.']
    longitude = df['Long.']
    age = df['Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]']
    
    # convert all data list to string
    index = index.astype(str)
    genetic_id = genetic_id.astype(str)
    countries_list = countries_list.astype(str)
    latitude = latitude.astype(str)
    longitude = longitude.astype(str)
    age = age.astype(str)
    
    sample_list = []
    for i in range(len(index)):
        # try block to catch any errors
        try:
            sample_list.append([index[i], genetic_id[i], countries_list[i], latitude[i], longitude[i], age[i]])
        except:
            pass

    # write the data to a new file
    with open(f'{path}/Ancient_samples.txt', 'w') as f:
        f.write('Index\tID\tCountry\tLatitude\tLongitude\tAge\n')
        for sample in sample_list:
            f.write('\t'.join(sample) + '\n')

#this function calculates the euclidean distance matrix from the plink output and writes it to a pickle file
def create_dist_matrix(df):
    # get the needed columns from the dataframe
    names = df['Genetic ID']
    lat = df['Lat.']
    long = df['Long.']
    admix = df.iloc[:,9:]
    
    # convert the latitude and longitude to float
    lat = lat.astype(float)
    long = long.astype(float)
    
    # covert data to numpy array
    names = np.array(names)
    lat = np.array(lat)
    long = np.array(long)
    admix = np.array(admix)
    
    #calculate the euclidean distance matrix
    dist = pdist(admix, metric='euclidean')
    dist_matrix = squareform(dist)
    dist_df = pd.DataFrame(dist_matrix, index=names, columns=names)
    
    # save the distance matrix to a pickle file
    dist_df.to_pickle("1_dist_matrix/eucl_dist.pkl")
    
# main function
if __name__ == "__main__":
    # check if the path to the excel file is given
    if len(sys.argv) != 2:
        print("Please provide the path to the excel file.")
        sys.exit()
    # check if the path is correct
    if not os.path.exists(sys.argv[1]):
        print("The path to the excel file is not correct.")
        sys.exit()
    # get the path to the excel file
    path = sys.argv[1]
    dict_path = os.path.dirname(path)
    
    #try to load the data into the dataframe
    try:
        df = pd.read_excel(path)
    except Exception as e:
        print("Could not read the excel file.")
        print(e)
        sys.exit()
        
    #filter the data to get rid of missing values for Latitude and Longitude
    df = filter_df(df)
    
    print("filtering ancient and modern samples...")
    write_sample_list(df, dict_path)
    
    print("calculating ibs distance matrix...")
    create_dist_matrix(df)
    
    