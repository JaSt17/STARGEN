
# Temporal Mapping of Genetic Barriers (TMGB)

## Overview

TMGB is a novel bioinformatics tool designed to visualize genetic distances between ancient DNA samples across different time periods and geographic locations. By integrating spatial and temporal dimensions, TMGB provides insights into genetic barriers, migration routes, and isolated populations, thereby contributing to the understanding of human population dynamics and historical gene flow.

## Prerequisites

Before running the tool, ensure you have the following software and libraries installed:

- Python 3.12.2
- conda 23.7.4

## Installation

1. Clone this repository to your local machine:

   ```bash
   git clone https://github.com/JaSt17/BINP17_TMGB.git
   cd BINP17_TMGB
   ```

2. Run the setup script to create a new conda environment and install all necessary packages.

   On Linux run:

   ```bash
   chmod +x setup.sh
   setup/./setup.sh
   ```

   On Windows run:

   ```powershell
   setup\.\setup.ps1
   ```

   If the installation of the conda environment and dependencies does not work using the setup file, please install the dependencies manually:

   - streamlit=1.32.0
   - streamlit-folium=0.20.0
   - folium=0.16.0
   - pandas=2.2.1
   - numpy=1.26.4
   - scipy=1.13.0
   - statsmodels=0.14.1
   - h3-py=3.7.6
   - haversine=2.8.1
   - branca=0.7.1
   - matplotlib=3.8.3
   - openpyxl=3.1.2

## Initial Run

Before running the tool for the first time, ensure you have the "aDNA_30GPs.xlsx" file in the `0_data/` directory. This file contains the ancient DNA samples necessary for the analysis.

Run the initial script to process the data and create the necessary distance matrix and other files:

On Linux run:

```bash
python scr/initial_run.py 0_data/aDNA_30GPs.xlsx
```

On Windows run:

```powershell
python scr\initial_run.py 0_data\aDNA_30GPs.xlsx
```

This script performs the initial data processing, which includes:

1. **Data Filtering**: Reads the Excel file and filters out rows with invalid or missing latitude and longitude values.
2. **Sample List Generation**: Extracts relevant columns from the filtered data and writes the ancient samples to a new text file.
3. **Distance Matrix Calculation**: Computes the Euclidean distance matrix from the admixture data and saves it as a pickle file.

### Output
- `Ancient_samples.txt`: A text file with the filtered ancient sample data.
- `1_dist_matrix/eucl_dist.pkl`: A pickle file containing the Euclidean distance matrix.

Ensure the path to the Excel file is correct when running the script to avoid errors.

## Pre-Execution Checklist

- Ensure the conda environment "tmbg" or your chosen environment is created and activated.
- Ensure your working directory contains the following files:
  - `0_data/Ancient_samples.txt`
  - `1_dist_matrix/eucl_dist.pkl`

## Running the Application

Start the interactive TMGB application with the following command:

On Linux run:

```bash
streamlit run scr/app.py
```

On Windows run:

```powershell
streamlit run scrpp.py
```

This command should automatically open a browser window with the TMGB application.

### Application Features

- **Home Screen**: Customize the number of time segments and the clarity of the hexagonal zones. Access additional details through the buttons on the home screen.

  ![Home Screen](./img/home_screen.png "Home Screen")

- **Main Screen**: Provides additional settings to adjust the visualization of genetic distances.

  ![Main Screen](./img/main_screen.png "Main Screen")

  - **Display Settings**:
    - **Which distances should be displayed**: Exclude distant hexagons below a set threshold.
    - **Isolated populations**: Define the scaled distance value to identify isolated populations and migration routes.
    - **Migration routes & isolated populations**: Draw possible migration routes and isolated hexagons on the map.
    - **Distance lines**: Display underlying distance lines for a detailed view.
    - **Default map window**: Set standard coordinates and zoom level for specific areas of interest.

Now, you can explore genetic distances over time using TMGB.

## Repository Structure

- **0_data**: Holds three xlsx files.
  - `aDNA_30GPs.xlsx`: Contains the combined samples from the AADR and additional samples gathered for this study.
  - `aDNA_30GPs_AADR.xlsx`: Contains only the AADR dataset samples.
  - `aDNA_30GPs_new_samples.xlsx`: Contains the new samples gathered for this study.

  Each of these files lists the Genetic ID, the Publication, the Date of the sample, the Locality, the Political entity, the coordinates, the SNPs hit on autosomal targets, and the 30 Admixture components.

- **scr**: Contains all the scripts necessary for the tool.
  - `app.py`: Holds the code for the Streamlit application. It allows the interactive usage of the tool.
  - `func.py`: Contains all the necessary functions to run the `app.py` script.
  - `initial_run.py`: Script to be run in the initial setup phase to create the distance matrix from the `aDNA_30GPs.xlsx` file.
  - `visualize.py`: Contains all functions used to visualize the hexagons, lines, and legend on the Folium map.

## Methodology and Validation

TMGB employs several advanced techniques for genetic data analysis, including:

- **Data Retrieval**: Uses datasets from the Allan Ancient DNA Resources (AADR) and other supplementary sources.
- **Euclidean Distance Computation**: Utilizes supervised ADMIXTURE analyses to calculate genetic distances.
- **Spatial Segmentation**: Implements the H3geo algorithm for hexagonal spatial segmentation.
- **Delaunay Triangulation**: Ensures robust spatial relationships between hexagons.
- **LOWESS Curve Fitting**: Scales genetic distances to geographic distances for each time bin.
- **Neighboring Hexagon Identification**: Identifies and visualizes genetic distances between hexagons.

TMGB has been validated by identifying known genetic barriers, such as the Bering Sea, and accurately detecting migration routes, such as those influenced by Norse maritime activities.

For detailed methodology and results, refer to the supplementary material and figures included in the manuscript. 

## References

A comprehensive list of references is included in the manuscript. Key references include foundational studies on genetic population structure, geographic patterns of genetic variation, and methodologies for spatial genetic analysis.

---

This README provides detailed instructions for running TMGB and replicating the results presented in the associated manuscript. For further information and updates, visit the [TMGB GitHub repository](https://github.com/JaSt17/BINP17_TMGB).
