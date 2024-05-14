# How to run the programm

## 1) Run the setup script

First we need to setup the conda environent and install all necessary dependencies.
To do so please run the setup script on the system of choice

On Linux run:

```{bash}
chmod +x setup.sh
./setup.sh
```

On Windows run:

```{powershell}
.\setup.ps1
```

## 2) Run the initial_run script

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
