if (!"already_run_00_wrapper" %in% ls()) {
  already_run_00_wrapper <- TRUE
  #######################################
  # Source functions
  fn_files_to_source <- list.files("functions/", full.names = TRUE)
  for (file_curr in fn_files_to_source) {
    source(file_curr)
  }
  control <- get_control_parameters_analysis()
  if (!"export_plots" %in% ls()) {
    export_plots <- FALSE
  }
  source("scripts/00_wrapper.R")
}


