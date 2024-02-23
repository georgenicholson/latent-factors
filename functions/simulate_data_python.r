simulate_data_python <- function(TASK_ID = 1, it = 1, 
                        sim_type = c("null", "non_null")[2]) {
    fn_files_to_source <- list.files("functions/", full.names = TRUE)
    for (file_curr in fn_files_to_source) {
       source(file_curr)
    }
    control <- get_control_parameters_simulations(parameter_set = TASK_ID)
    sim_curr <- simulate_data(control = control, 
                              seed = it + (TASK_ID - 1) * control$n_its, 
                              zero_ATE = ifelse(sim_type == "null", TRUE, FALSE))
   return(sim_curr)
}




#   placebo_label <- "Placebo"
#   estimate_type <- c("treatment_means", "contrasts_vs_placebo")[2]
#   var_names = paste0("Y", 1:control$n_var)
#   n_fac <- length(control$endpoint_proportions)
#   model_names <- c("fac_0_longit_0", "fac_0_longit_1", "fac_1_longit_0", "fac_1_longit_1")
#   restab <- data.frame(it = 1:control$n_its, est = NA, low = NA, upp = NA, se = NA, pval = NA)
#   resl <- list()
#   for (sim_type in c("non_null", "null")) {
#     resl[[sim_type]] <- list()
#     for (model_name in model_names) {
#     resl[[sim_type]][[model_name]] <- list()
#       resl[[sim_type]][[model_name]] <- list()
#       resl[[sim_type]][[model_name]]$y <- restab
#       for (fac_curr in 1:n_fac) {
#         resl[[sim_type]][[model_name]][[paste0("z", fac_curr)]] <- restab
#       }
#     }
#   }
#   for (it in 1:control$n_its) {
#     if (it %% 25 == 0) {
#       print(it)
#     }
  
#     sim_type <- "non_null"
#     for (sim_type in c("null", "non_null")) {
#       ##########################################
#       # Simulate data 
#       sim_curr <- simulate_data(control = control, 
#                               seed = it + (TASK_ID - 1) * control$n_its, 
#                               zero_ATE = ifelse(sim_type == "null", TRUE, FALSE))
    
    
# }