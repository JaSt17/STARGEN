#!/bin/bash

# Name of the environment
ENV_NAME="myenv"

# List of packages 
CONDA_PACKAGES=(
    "python=3.12.2"
    "streamlit=1.32.0"
    "streamlit-folium=0.19.0"
    "folium=0.16.0"
    "pandas=2.2.1"
    "numpy=1.26.4"
    "scipy=1.12.0"
    "h3-py=3.7.7"
    "matplotlib=3.8.0"
    "openpyxl=3.1.2"
)

# Create the Conda environment with the correct Python version
echo "Creating environment '$ENV_NAME' with Python version ${CONDA_PACKAGES[0]}"
conda create --name "$ENV_NAME" "${CONDA_PACKAGES[0]}" -y -c conda-forge

# Activate the Conda environment
echo "Activating environment '$ENV_NAME'"
conda activate "$ENV_NAME" 

# Install each package using Conda, starting from the second item
for (( i=1; i<${#CONDA_PACKAGES[@]}; i++ ))
do
    echo "Installing ${CONDA_PACKAGES[i]} using Conda..."
    conda install -n "$ENV_NAME" "${CONDA_PACKAGES[i]}" -y -c conda-forge
done

echo "All packages have been installed successfully."
