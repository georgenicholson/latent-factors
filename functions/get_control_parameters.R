#########################################################################
# Function to get plotting parameters
get_control_parameters <- function(parameter_set = 1) {
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

  
  control$parameter_change_list <- list("Default parameters" = list(),
                                        "Increased magnitude of (negative) effect size $\\mu_{\\text{active,16}}$" = list(size_of_final_week_effect_active = -1.4),
                                        "Decreased effect size" = list(size_of_final_week_effect_active = -1.1),
                                        "Decreased $\\rhoARNok$" = list(cor_ar1_lag = .25),
                                        "Increased $\\rhoARNok$" = list(cor_ar1_lag = .95),
                                        # Note the first element in each of the below lists represents the changed parameter,  
                                        # so the ordering below is important!!
                                        "Increased $\\sigmaARNok$" = list(theta_resid_sd = 1, theta_RI_sd = .5, y_resid_sd = .5),
                                        "Increased $\\sigmaRINok$" = list(theta_RI_sd = 1, theta_resid_sd = .5, y_resid_sd = .5),
                                        "Increased $\\sigmaY$" = list(y_resid_sd = 1, theta_resid_sd = .5, theta_RI_sd = .5))
  
  control$parameter_change_list
  change_curr <- control$parameter_change_list[[parameter_set]]
  for (nam in names(change_curr)) {
    control[[nam]] <- change_curr[[nam]]
  }
  control$parameter_set <- parameter_set
  control$n_parameter_sets <- length(control$parameter_change_list)
  return(control)
}




