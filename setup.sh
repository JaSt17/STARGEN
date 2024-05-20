#!/bin/bash

# Name of the environment
ENV_NAME="myenv"

# Required Conda version
REQUIRED_CONDA_VERSION="23.7.4"

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

# Function to check if Conda exists and its version
check_conda() {
    if ! command -v conda &> /dev/null
    then
        echo "Conda not found. Please install Conda version $REQUIRED_CONDA_VERSION."
        exit 1
    fi

    CONDA_VERSION=$(conda --version | awk '{print $2}')
    if [ "$CONDA_VERSION" != "$REQUIRED_CONDA_VERSION" ]; then
        echo "Warning: Conda version $CONDA_VERSION found. Please use Conda version $REQUIRED_CONDA_VERSION to avoid potential errors."
    fi
}

# Check Conda installation and version
check_conda

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
