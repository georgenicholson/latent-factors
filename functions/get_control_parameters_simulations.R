#########################################################################
# Function to get plotting parameters
get_control_parameters_simulations <- function(parameter_set = 1) {
  control <- list(

    # Default parameters for simulation (to be varied)
    n_subj = 200,
    cor_ar1_lag = .75,
    theta_RI_sd = 0.5,
    theta_resid_sd = 0.5,
    y_resid_sd = 0.5, 
    size_of_final_week_effect_active = -1.2,
    missing_proportion = 0, 
    set_baseline_arm_effects_equal = T,
    use_empirical_loadings = T, 
    
    # Default parameters for simulation (fixed)
    endpoint_proportions = c(5, 5, 2) / 12,
    n_its = 1000, # 1000 its takes 56 mins
    n_its_benchmark = 1000, # Fewer data sets analysed in benchmarking as external methods (RNN, VAR etc.) are more computationally intensive
    n_var = 12,
    n_fac = 3,
    time_of_half_max_effect = 2,
    size_of_final_week_effect_placebo = -1,
    weeks_measure = c(0, 1, 2, 4, 8, 12, 16),
    time_pt_eval = 16,
    ci_alpha = 0.05,
    pval_thresh = 0.05,
    y_ind_use = 1,
    TRT_prob = 0.5
  )

  
  control$parameter_change_list <- list("$\\Theta_1$: Default parameters" = list(),
                                        "$\\Theta_2$: Increased magnitude of (negative) effect size $\\mu_{\\text{active,16}}$" = list(size_of_final_week_effect_active = -1.4),
                                        "Decreased effect size" = list(size_of_final_week_effect_active = -1.1),
                                        "$\\Theta_4$: Decreased $\\rhoARNok$" = list(cor_ar1_lag = .25),
                                        "$\\Theta_3$: Increased $\\rhoARNok$" = list(cor_ar1_lag = .95),
                                        # Note the first element in each of the below lists represents the changed parameter,  
                                        # so the ordering below is important!!
                                        "$\\Theta_5$: Increased $\\sigmaARNok$" = list(theta_resid_sd = 1, theta_RI_sd = .5, y_resid_sd = .5),
                                        "$\\Theta_6$: Increased $\\sigmaRINok$" = list(theta_RI_sd = 1, theta_resid_sd = .5, y_resid_sd = .5),
                                        "$\\Theta_7$: Increased $\\sigmaY$" = list(y_resid_sd = 1, theta_resid_sd = .5, theta_RI_sd = .5))
  
  control$ts_methods <- c("car1", "arma_1_0", "arma_0_1", "arma_1_1", "arma_2_0", "arma_0_2", "arma_2_1", "arma_1_2", "exp_smooth")
  inds_keep <- c(1:8, 11:17)
  control$model_names <- c("fac_0_longit_0", 
                           "fac_1_longit_0", 
                           paste0("fac_0_longit_1_", control$ts_methods), 
                           "fac_1_longit_1",
                           "LM4_linear",
                           "VAR",
                           "RNN",
                           "LSTM",
                           "GRU")[inds_keep]
  control$input <- c(1, 2, rep(3, length(control$ts_methods)), rep(4, 6))[inds_keep]
  control$methods <- c("\\lm{1}", 
                       "\\lm{2}",
                       "\\lm{3}", 
                       "ARMA(1,0)", 
                       "ARMA(0,1)", 
                       "ARMA(1,1)", 
                       "ARMA(2,0)", 
                       "ARMA(0,2)", 
                       "ARMA(2,1)", 
                       "ARMA(1,2)",
                       "ES",
                       "\\lm{4}",
                       "\\lmlinear{4}",
                       "VAR",
                       "RNN",
                       "LSTM",
                       "GRU")[inds_keep]
  
  control$complexity <- c("\\lm{1}" = "$O(PN)$",
                          "\\lm{2}" = "$O(PK^2 + PNT + KN)$",
                          "CAR(1) (\\lm{3})" = "$O(PNT^3)$",
                          "ARMA(1,0)" = "$O(PNT^3)$",
                          "ARMA(0,1)" = "$O(PNT^3)$",
                          "ARMA(1,1)" = "$O(PNT^3)$",
                          "ARMA(2,0)" = "$O(PNT^3)$",
                          "ARMA(0,2)" = "$O(PNT^3)$",
                          "ARMA(2,1)" = "$O(PNT^3)$",
                          "ARMA(1,2)" = "$O(PNT^3)$",
                          "ES" = "$O(PNT^3)$",
                          "\\lm{4}" = "$O(PK^2 + PNT + KNT^3)$",
                          "\\lmlinear{4}" = "$O(PK^2 + PNT + KNT^3)$",
                          "VAR" = "$O(N P^3 + N P^2 T)$", 
                          "RNN" = "$O(NTE(PH + H^2))$", 
                          "LSTM" = "$O(NTE(PH + H^2))$", 
                          "GRU" = "$O(NTE(PH + H^2))$")[inds_keep]
  control$parameter_change_list
  change_curr <- control$parameter_change_list[[parameter_set]]
  for (nam in names(change_curr)) {
    control[[nam]] <- change_curr[[nam]]
  }
  control$parameter_set <- parameter_set
  control$n_parameter_sets <- length(control$parameter_change_list)
  
  control$sim_results_file <- file.path("tables", "sim_results.txt")
  control$default_sim_table_file <- file.path("tables", "default_sim_parameter_table.txt")
  control$emax_example_plot_file <- file.path("plots", "example_emax_curves.pdf")
  # control$benchmarking_results_file <- file.path("tables", "benchmarking_results.txt")
  control$benchmarking_runtime_file <- file.path("tables", "benchmarking_runtime_results.txt")
  control$numbers_out_dir <- "text_numbers/"
  control$benchmarking_table_split <- list(c(1, 2), c(5, 4), 6:8)
  control$benchmarking_results_files <- c()
  for(i in 1:length(control$benchmarking_table_split)) {
    control$benchmarking_results_files[i] <- file.path("tables", paste0("benchmarking_results_par_settings_", 
                                                                        paste(control$benchmarking_table_split[[i]], collapse = "-"), 
                                                                        ".txt"))
  }
  
  return(control)
}




