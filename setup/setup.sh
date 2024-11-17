#!/bin/bash

# Name of the environment
ENV_NAME="stargen"

# List of packages 
CONDA_PACKAGES=(
    "python=3.12.2"
    "streamlit=1.32.0"
    "streamlit-folium=0.20.0"
    "folium=0.16.0"
    "pandas=2.2.1"
    "numpy=1.26.4"
    "scipy=1.13.0"
    "statsmodels=0.14.1"
    "h3-py=3.7.6"
    "haversine=2.8.1"
    "branca=0.7.1"
    "matplotlib=3.8.3"
    "openpyxl=3.1.2"
)

# Function to check if Conda exists
check_conda() {
    if ! command -v conda &> /dev/null
    then
        echo "Conda could not be found. Please install Conda first."
        exit 1
    fi
}

# Check Conda installation
check_conda

# Create the Conda environment with the correct Python version
echo "Creating environment '$ENV_NAME' with Python version ${CONDA_PACKAGES[0]}"
conda create --name "$ENV_NAME" "${CONDA_PACKAGES[0]}" -y -c conda-forge

# Activate the Conda environment
echo "Activating environment '$ENV_NAME'"
source $(conda info --base)/etc/profile.d/conda.sh
conda activate "$ENV_NAME"

# Install each package using conda
for package in "${CONDA_PACKAGES[@]:1}"
do
    echo "Installing ${package} using conda..."
    conda install "$package" -y -c conda-forge
done

echo "All packages have been installed successfully."

# Activate the Conda environment
conda activate "$ENV_NAME"

echo "Environment '$ENV_NAME' is ready to use."
