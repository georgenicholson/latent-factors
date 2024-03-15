#######################################
# Load local R package
devtools::load_all("mlfa_package")

#######################################
# Source functions
fn_files_to_source <- list.files("functions/", full.names = TRUE)
for (file_curr in fn_files_to_source) {
  source(file_curr)
}


sim_type <- "non_null"
n_its_for_python <- 1000
for (TASK_ID in 1:8) {
  print(TASK_ID)
  control <- get_control_parameters_simulations(parameter_set = TASK_ID)
  dir_curr <- file.path("sims_out", "CSVs_sim_", TASK_ID)
  dir.create(dir_curr, showWarnings = FALSE, recursive = TRUE)
  for (it in 1:n_its_for_python) {
    print(it)
    sim_out <- simulate_data(control = control, 
                             seed = it + (TASK_ID - 1) * control$n_its, 
                             zero_ATE = ifelse(sim_type == "null", TRUE, FALSE))
    write.csv(sim_out$Y, 
              file = file.path(dir_curr, paste0("Y_", it, ".csv")), 
              row.names = FALSE)
  }
}  


