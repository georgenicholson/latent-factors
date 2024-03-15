#######################################
# Load local R package
devtools::load_all("mlfa_package")

#######################################
# Source functions
fn_files_to_source <- list.files("functions/", full.names = TRUE)
for (file_curr in fn_files_to_source) {
  source(file_curr)
}

#######################################
# Folder structure
DIR <- get_DIR()

computation_type <- c("simulations", "data_analysis")[1]

#######################################
# Simulations
if (computation_type == "simulations") {
  source("scripts/s01_emax_curves_demo_plot.R")
  source("scripts/s02_simulation_study.R")
  source("scripts/s03_output_latex_tables.R")
  source("scripts/s04_benchmarking.R")
  source("scripts/s05_output_benchmarking_latex_tables.R")
  source("scripts/s06_export_results.R")
}

#######################################
# Data analysis (this requires application for data access and import)
if (computation_type == "data_analysis") {
  
  force_load_data <- FALSE
  force_preproc_data <- FALSE
  force_analyse_data <- FALSE
  
  #######################################
  # Load data
  if (!"d_step_1.RDS" %in% list.files("preproc_data") | force_load_data) {
    source("scripts/01_load_PsA_RA_newRelease.R")
  } else {
    d_step1 <- readRDS(file = "../preproc_data/d_step_1.RDS")
    data_vars_all <- readRDS(file = "../preproc_data/data_vars_all.RDS")
    d_subject_00 <- readRDS(file = "../preproc_data/d_subject_00.RDS")
  }
  
  #######################################
  # Preprocess data
  if (!"d_time.RDS" %in% list.files("preproc_data") | force_preproc_data) {
    source("scripts/02_preprocess_PsA_RA.R")
  } else {
    d_time <- readRDS(file = "../preproc_data/d_time.RDS")
    d_subject <- readRDS(file = "../preproc_data/d_subject.RDS")
    nms_params_v2 <- readRDS(file = "../preproc_data/nms_params_v2")
  }
  endpts <- nms_params_v2$PARAM
  n_endpt <- length(endpts)
  
  #######################################
  # Analyse the data
  if (!"analysis_results_list.RDS" %in% list.files("output") | force_analyse_data) {
    source("scripts/03_analysis.R")
  } else {
    analysis_results_list <- readRDS(file = file.path("output", "analysis_results_list.RDS"))
  }
  
  #######################################
  # Generate plots
  source("scripts/04z_run_paper_plot_scripts.R")
}

