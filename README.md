# Temporal Mapping of Genetic Barriers (TMGB)

## Overview

This README file is a step to step intruction on how to run the tool and replicate the results presented in the Manuscript. The tool is designed to perform visualise gentic distance between Ancient DNA Samples in different areas of time.

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

2. Run the setup scipt to create a new conda env and install all nedded packages

  First we need to setup the conda environent and install all necessary dependencies.
  To do so please run the setup script on the system of choice

  On Linux run:

   ```{bash}
    chmod +x setup.sh
    setup/./setup.sh
   ```

  On Windows run:

   ```{powershell}
    setup\.\setup.ps1
   ```

  If the installation of the conda envirometn and dependencies should not work using the setup file please install the dependecies needed (listed in `setup/requirements.txt`) manually.
3. Run the initial_run script

Before running the tool for the first time make sure to you have the "aDNA_30GPs.xlsx" file in the 0_data/ directory.
Than run the the initial run.

On Linux run:

```{bash}
python scr/initial_run.py 0_data/aDNA_30GPs.xlsx 
```

On Windows run:

```{powershell}
python scr\initial_run.py 0_data\aDNA_30GPs.xlsx 
```

## 3) Before running

- make sure that you created and activated the conda environment "myenv"

- make sure that your working directory contains the following files:
  - 0_data/Ancient_samples.txt
  - 1_dist_matrix/eucl_dist.pkl

## 4) Run the application

Now the interactive application can be started with the following bash command.

On Linux run:

```{bash}
streamlit run scr/app.py
```

On Windows run:

```{powershell}
streamlit run scr\app.py
```
