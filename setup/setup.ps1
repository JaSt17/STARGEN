# Name of the environment
$ENV_NAME = "stargen"

# List of Conda packages
$CONDA_PACKAGES = @(
    "python=3.12.2",
    "streamlit=1.32.0",
    "streamlit-folium=0.20.0",
    "folium=0.16.0",
    "pandas=2.2.1",
    "numpy=1.26.4",
    "scipy=1.13.0",
    "statsmodels=0.14.1",
    "h3-py=3.7.6",
    "haversine=2.8.1",
    "branca=0.7.1",
    "matplotlib=3.8.3",
    "openpyxl=3.1.2"
)

# Function to check if Conda exists and its version
function Check-Conda {
    if (-not (Get-Command conda -ErrorAction SilentlyContinue)) {
        Write-Host "Conda not found. Please install Conda first."
        exit 1
    }
}

# Check Conda installation and version
Check-Conda

# Create the Conda environment with Python version
Write-Host "Creating environment '$ENV_NAME' with Python version $($CONDA_PACKAGES[0])"
conda create --name $ENV_NAME $($CONDA_PACKAGES[0]) -y -c conda-forge

# Activate the Conda environment
Write-Host "Activating environment '$ENV_NAME'"
conda activate $ENV_NAME

# Ensure the conda environment activation is available in the script
& $env:CONDA_PREFIX\Scripts\activate.ps1

# Install each package using Conda, starting from the second item
foreach ($package in $CONDA_PACKAGES[1..($CONDA_PACKAGES.Length - 1)]) {
    Write-Host "Installing $package using Conda..."
    conda install -n $ENV_NAME $package -y -c conda-forge
}

Write-Host "All packages have been installed successfully."

# activate the environment
conda activate $ENV_NAME

Write-Host "Environment '$ENV_NAME' is ready to use."
