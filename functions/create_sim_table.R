create_sim_table <- function(n_subj = c(50, 100, 200),
                             theta_resid_sd = c(1, 0),
                             cor_ar1_lag = c(.99, .75, 0),
                             theta_RI_sd = 1,
                             y_resid_sd = c(.5, 1, 2),
                             endpoint_proportions = list(c(0.6, 0.3, 0.1),
                                                         c(1, 1, 1, 1) / 4),
                             missing_proportion = c(0, .25),
                             beta_Emax = c(0, .1, .25, .5),
                             set_baseline_arm_effects_equal = c(TRUE, FALSE),
                             TRT_prob = 0.5,
                             E0_null = 0,
                             E50_null = 7,
                             Emax_null = 3,
                             beta_E0 = 0,
                             beta_E50 = 0,
                             n_var = 12,
                             weeks_measure = c(0, 1, 2, 4, 8, 12, 16),
                             time_pt_eval = "16") {



}

