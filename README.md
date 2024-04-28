# A framework for longitudinal latent factor modelling of treatment response in clinical trials with applications to Psoriatic Arthritis and Rheumatoid Arthritis

This repository contains the R scripts needed to reproduce the simulated results reported 
in the manuscript 'A framework for longitudinal latent factor modelling of treatment response in clinical trials with applications to Psoriatic Arthritis and Rheumatoid Arthritis'. 

## Installation

To run these R scripts, you will need R version 3.6.3 or later, widely available on 
Unix-like, Windows and Mac families of operating systems. If you have R version 4.0.0
or above on Windows, you will also need to install 
[Rtools](https://cran.r-project.org/bin/windows/Rtools/). To get started,
first [clone](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository)
this repository onto your local machine. Next, open an R console, and install the 
[renv](https://rstudio.github.io/renv/index.html) R package if you don't have it 
already (e.g. via `install.packages("renv")`). Then, run the following code block to install the required packages for the scripts, changing `path_to_dir` to the path of your local version of this repository (note that you need to change it twice). 

```
path_to_dir <- "path/to/latent-factors"
setwd(path_to_dir)
renv::activate(path_to_dir)
path_to_dir <- "path/to/latent-factors"
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
renv::restore(path_to_dir)
```

In the final step running `renv::restore(path_to_dir)`, make sure to select option 2 ("2: Do not activate the project and use the current library paths.").

- - -

To run the Python scripts, you need to install a virtual or conda environment with the dependencies as listed in `requirements.txt`.

Further, you need to download the simulation datasets [here](https://drive.google.com/file/d/1YzcPeoYiNyJLaUKkWHyouDLr_2omVw6x/view?usp=sharing). Unzip the folder, and place it in `data/`.


## Scripts to reproduce simulation results in manuscript

The subdirectory `scripts` contains all the code needed to reproduce the results. 
This can be conveniently done by sourcing the wrapper function as follows:

```
source("scripts/00_wrapper.R")
```

The above script will first source `scripts/s01_emax_curves_demo_plot.R`, creating a plot of example Emax curves in `plots/example_emax_curves.pdf`.

Next the script will source  `scripts/s02_simulation_study.R`, which will run the simulation study.

Finally the script sources  `scripts/s03_output_latex_tables.R`, which collects results and outputs tables in Latex format in `tables/default_sim_parameter_table.txt`.

- - -

Before we begin, at the top of `p09_rnn.py`, `p09_eval_rnn.py`, `p09_attach_sweep.py` and `p09_start_sweep.py`, change the `project` and `entity` variable to your `wandb` settings.
Further, at the top of `p09_rnn.py`, change `repo_path` to the directory where the `latent_factors` repo is located.

To reproduce the RNN experiments, run a sweep.

```python
python3 p09_start_sweep.py
```

This uses the simulation data in `data/sims_out/`.

You may also attach more runs to the sweep using the `p09_attach_sweep.py` script. 

Then, change the variable `sweep_id` with the corresponding `wandb` sweep ID, and run the script `p09_eval_rnn.py` to compute the results.
