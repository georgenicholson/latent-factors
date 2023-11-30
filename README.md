# A framework for longitudinal latent factor modelling of treatment response in clinical trials with applications to Psoriatic Arthritis and Rheumatoid Arthritis

This repository contains the R scripts needed to reproduce the simulated results reported 
in the manuscript 'A framework for longitudinal latent factor modelling of treatment response in clinical trials with applications to Psoriatic Arthritis and Rheumatoid Arthritis'. 

## Installation

To run these scripts, you will need R version 3.6.3 or later, widely available on 
Unix-like, Windows and Mac families of operating systems. If you have R version 4.0.0
or above on Windows, you will also need to install 
[Rtools](https://cran.r-project.org/bin/windows/Rtools/). To get started,
first [clone](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository)
this repository onto your local machine. Next, open an R console, and install the 
[renv](https://rstudio.github.io/renv/index.html) R package if you don't have it 
already (e.g. via `install.packages("renv")`). Then, run the following code block to install the required packages for the scripts,  
changing `path_to_dir` to the path of your local version of this repository (note that you need to change it twice). 

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


## Scripts to reproduce simulation results in manuscript

The subdirectory `scripts` contains all the code needed to reproduce the results. 
This can be conveniently done by sourcing the wrapper function as follows:

source("scripts/00_wrapper.R")
