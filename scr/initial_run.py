# This script is responsible for the initial run of the pipeline.
# It reads the excel file with the samples, filters the data, and writes the samples to a new file.
# It also calculates the euclidean distance matrix from the plink output and writes it to a pickle file.
# The script is called with the path to the excel file as an argument.

#-------------------------------------------------------------------------------------
# NOTE: If you want to use this script with a different Excel file,                  
# you need to adjust the columns list underneath                
columns = ['Genetic ID',
                'Lat.',
                'Long.', 
                'Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]']
# also change the index where the the distance values begins
index = 9 # Assuming data starts from the 10th column
#-------------------------------------------------------------------------------------

import pickle
import os
import sys
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform

def filter_df(df):
    """
    Filters the DataFrame for errors in the latitude and longitude columns.

    Parameters:
    - df (pd.DataFrame): The input DataFrame containing 'Lat.' and 'Long.' columns.

    Returns:
    - pd.DataFrame: The filtered DataFrame with valid latitude and longitude values.
    """
    # Filter out samples with missing or invalid latitude and longitude values
    filtered_df = df[(df['Lat.'] != '..') & (df['Long.'] != '..')]
    
    # Drop rows with NaN values in 'Lat.' and 'Long.' columns
    filtered_df = filtered_df.dropna(subset=['Lat.', 'Long.'])
    
    return filtered_df


def write_sample_list(df, path, columns):
    """
    Writes the samples to a new file.

    Parameters:
    - df (pd.DataFrame): The input DataFrame containing the sample data (needs to include 'Genetic ID', 'Lat.', 'Long.', and 'Date mean in BP' columns)
    - path (str): The path to the directory where the file will be saved.

    Returns:
    - None
    """
    # Extract the required columns and convert them to string
    # try to find the columns in the dataframe and return an error if no match:
    try:
        data = df[columns].astype(str)
    except KeyError as e:
        print("ERROR:\nCould not find all necessary columns in your excel file.\nPlease change the names of your columns according to your provided excel file in the initial_run.py script.")
        sys.exit(1)

    # Prepare the sample list
    sample_list = data.values.tolist()

    # Write the data to a new file
    file_path = f'{path}/Ancient_samples.txt'
    with open(file_path, 'w') as f:
        f.write('ID\tLatitude\tLongitude\tAge\n')
        for sample in sample_list:
            f.write('\t'.join(sample) + '\n')


def create_dist_matrix(df, columns, index):
    """
    Calculates the Euclidean distance matrix from the samples and writes it to a pickle file.

    Parameters:
    - df (pd.DataFrame): The input DataFrame containing 'Genetic ID', 'Lat.', 'Long.', and admixture columns.

    Returns:
    - None
    """
    # Extract the needed columns from the DataFrame
    names = df[columns[0]]
    lat = df[columns[1]]
    long = df[columns[2]]
    try:
        admix = df.iloc[:, index:]  # Assuming data starts from the 10th column
        admix.astype(float)
    except:
        print("There is something wrong with the index you set for the distance values.\nPlease check your provided index in the intial_run.py file.")
        sys.exit(1)
    # check if the distance values dataframe is empty and if so stop the process and let the user know.
    if admix.empty:
        print("Could not find distance values.\nPlease check your provided index in the intial_run.py file.")
        sys.exit(1)
    # tell the user the size of their distance vector so they can adjsut it if it is not correct.
    print(f"Your distance vector holds {admix.shape[1]} values.\nIf that is not correct, please check your provided index in the intial_run.py file. ")
    # Convert latitude and longitude to float
    lat = lat.astype(float)
    long = long.astype(float)
    
    # Convert data to numpy arrays
    names = np.array(names)
    lat = np.array(lat)
    long = np.array(long)
    admix = np.array(admix)
    
    # Calculate the Euclidean distance matrix
    dist = pdist(admix, metric='euclidean')
    dist_matrix = squareform(dist)
    dist_df = pd.DataFrame(dist_matrix, index=names, columns=names)
    
    # Create output directory if it doesn't exist
    os.makedirs("1_dist_matrix", exist_ok=True)
    
    # Save the distance matrix to a pickle file
    dist_df.to_pickle("1_dist_matrix/eucl_dist.pkl")


def main():
    """
    Main function to process the Excel file, filter data, write sample list, and calculate distance matrix.
    
    Usage:
    python script_name.py <path_to_excel_file>
    """
    # Check if the path to the Excel file is given
    if len(sys.argv) != 2:
        print("Please provide the path to the Excel file.")
        sys.exit(1)
    
    # Get the path to the Excel file
    path = sys.argv[1]
    
    # Check if the path is correct
    if not os.path.exists(path):
        print("The path to the Excel file is not correct.")
        sys.exit(1)
    
    dict_path = os.path.dirname(path)
    
    # Try to load the data into the DataFrame
    try:
        df = pd.read_excel(path)
    except Exception as e:
        print("Could not read the Excel file.")
        print(e)
        sys.exit(1)
        
    # Filter the data to get rid of missing values for Latitude and Longitude
    df = filter_df(df)
    
    # Process the filtered data
    print("Filtering ancient and modern samples...")
    write_sample_list(df, dict_path, columns)
    
    print("Calculating euclidean distance matrix...")
    create_dist_matrix(df, columns, index)
    
    print("Finished succesfully!")

if __name__ == "__main__":
    main()