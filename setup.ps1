# Name of the environment
$ENV_NAME = "myenv"

# Required Conda version
$REQUIRED_CONDA_VERSION = "23.7.4"

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
        Write-Host "Conda not found. Please install Conda version $REQUIRED_CONDA_VERSION."
        exit 1
    }

    $condaVersionOutput = conda --version
    $condaVersion = $condaVersionOutput -replace "conda ", ""
    if ($condaVersion -ne $REQUIRED_CONDA_VERSION) {
        Write-Host "Warning: Conda version $condaVersion found. Please use Conda version $REQUIRED_CONDA_VERSION to avoid potential errors."
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

# Install each package using Conda, starting from the second item
foreach ($package in $CONDA_PACKAGES[1..($CONDA_PACKAGES.Length - 1)]) {
    Write-Host "Installing $package using Conda..."
    conda install -n $ENV_NAME $package -y -c conda-forge
}

Write-Host "All packages have been installed successfully."
