# Name of the environment
$ENV_NAME = "myenv"

# List of Conda packages
$CONDA_PACKAGES = @(
    "python=3.12.2",
    "streamlit=1.32.0",
    "streamlit-folium=0.19.0",
    "folium=0.16.0",
    "pandas=2.2.1",
    "numpy=1.26.4",
    "scipy=1.12.0",
    "h3-py=3.7.7",
    "matplotlib=3.8.0",
    "openpyxl=3.1.2", 
    "pykrige=1.7.1"
)

# Create the Conda environment with Python version
Write-Host "Creating environment '$ENV_NAME' with Python version $($CONDA_PACKAGES[0])"
conda create --name $ENV_NAME $($CONDA_PACKAGES[0]) -y -c conda-forge

# Activate the Conda environment
Write-Host "Activating environment '$ENV_NAME'"
conda activate $ENV_NAME

# Install each package using Conda, starting from the second item
foreach ($package in $CONDA_PACKAGES[1..$CONDA_PACKAGES.Length]) {
    Write-Host "Installing $package using Conda..."
    conda install -n $ENV_NAME $package -y -c conda-forge
}

Write-Host "All packages have been installed successfully."
