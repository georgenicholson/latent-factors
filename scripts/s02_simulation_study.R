.libPaths("/gpfs3/well/lindgren-ukbb/users/qgj547/R/4.2/skylake")

#######################################
# Load local R package
devtools::load_all("mlfa_package")

#######################################
# Source functions
fn_files_to_source <- list.files("functions/", full.names = TRUE)
for (file_curr in fn_files_to_source) {
  source(file_curr)
}

run_on_a_server <- FALSE
if (run_on_a_server) {
  source("functions/process_command_line_args.R")
  args <- process_command_line_args(args_default = list(TASK_ID = 1))
  TASK_ID_ALL <- args$TASK_ID
} else {
  TASK_ID_ALL <- 1:8
}

dir.create("output", showWarnings = FALSE)
for (TASK_ID in TASK_ID_ALL) {
  control <- get_control_parameters_simulations(parameter_set = TASK_ID)
  placebo_label <- "Placebo"
  estimate_type <- c("treatment_means", "contrasts_vs_placebo")[2]
  var_names = paste0("Y", 1:control$n_var)
  n_fac <- length(control$endpoint_proportions)
  model_names <- c("fac_0_longit_0", "fac_0_longit_1", "fac_1_longit_0", "fac_1_longit_1")
  restab <- data.frame(it = 1:control$n_its, est = NA, low = NA, upp = NA, se = NA, pval = NA)
  resl <- list()
  for (sim_type in c("non_null", "null")) {
    resl[[sim_type]] <- list()
    for (model_name in model_names) {
    resl[[sim_type]][[model_name]] <- list()
      resl[[sim_type]][[model_name]] <- list()
      resl[[sim_type]][[model_name]]$y <- restab
      for (fac_curr in 1:n_fac) {
        resl[[sim_type]][[model_name]][[paste0("z", fac_curr)]] <- restab
      }
    }
  }
  for (it in 1:control$n_its) {
    if (it %% 25 == 0) {
      print(it)
    }
    for (sim_type in c("null", "non_null")) {
      ##########################################
      # Simulate data 
      sim_curr <- simulate_data(control = control, 
                              seed = it + (TASK_ID - 1) * control$n_its, 
                              zero_ATE = ifelse(sim_type == "null", TRUE, FALSE))
    
    
      if (sim_type == "non_null") {
        sim_non_null <- sim_curr
      } else {
        sim_null <- sim_curr
      }
      data_in <- as_tibble(sim_curr$Y)
      sub_data <- list(time = as_tibble(sim_curr$Y),
                       subject = as_tibble(sim_curr$dbase))
  
      ##########################################
      # Stage 1: estimate sparse loadings A and scores z
      mfd_object <- mlfa::mfd(data = sub_data$time, 
                  vars = var_names, 
                  loadings.n_PCs = n_fac, 
                  time = "Time", 
                  subject = "ID",
                  normalize.type = c("per_meas_standard", "standardize", "per_center_baseline_scale")[1],
                  Time_base = 0, 
                  standardize_A_type = c("none", "L2", "largest_abs")[3],
                  loadings.sparse_abs_thres = 0,
                  loadings.rotation = c("varimax", "promax")[1],
                  scores.type = 'a_posteriori', 
                  trace = FALSE)
  
       
      ##########################################
      # Assemble data for Stage 2
      arm_all <- unique(sub_data$subject$TRT)
      n_arm <- length(arm_all)
      results_df <- as.data.frame(mfd_object$data_out)
      y_nam <- paste0("Y", control$y_ind_use)
      results_df$arm <- sub_data$subject$TRT[match(results_df$ID, sub_data$subject$ID)]
      results_df$arm <- relevel(factor(results_df$arm), ref = "Placebo")
      results_df$time_fac <- factor(results_df$Time)
      results_df$id_fac <- factor(results_df$ID)
      nt <- length(control$weeks_measure)
      m_coef_mat<-u_coef_mat<-l_coef_mat <- matrix(NA, nt, n_arm, dimnames = list(control$weeks_measure, arm_all))
      results_list <- list()
      results_df$y <- sub_data$time %>% pull(paste0("Y", control$y_ind_use))
      for (time_curr in control$weeks_measure) {
        results_df[, paste0("time_", time_curr, "_binary")] <- ifelse(results_df$Time == time_curr, 1, 0)
      }
      non_baseline_form <- paste(paste0("time_", setdiff(control$weeks_measure, 0), "_binary"), collapse = " + ")
      include_baseline_form <- paste(paste0("time_", control$weeks_measure, "_binary"), collapse = " + ")
      if (control$set_baseline_arm_effects_equal) {
        baseline_form_curr <- non_baseline_form
      } else {
        baseline_form_curr <- include_baseline_form
      }
      results_df$active_binary <- ifelse(results_df$arm == "Active", 1, 0)
      
      ##########################################################
      # Single time point analysis of single measurement
      active_y_single_t <- results_df[results_df$arm == "Active" & results_df$Time == control$time_pt_eval, "y"]
      placebo_y_single_t <- results_df[results_df$arm == "Placebo" & results_df$Time == control$time_pt_eval, "y"]
      t_test_curr <- t.test(x = active_y_single_t, y = placebo_y_single_t, var.equal = TRUE)
      resl[[sim_type]]$fac_0_longit_0$y$est[it] <- t_test_curr$estimate[1] - t_test_curr$estimate[2]
      resl[[sim_type]]$fac_0_longit_0$y$low[it] <- t_test_curr$conf.int[1]
      resl[[sim_type]]$fac_0_longit_0$y$upp[it] <- t_test_curr$conf.int[2]
      resl[[sim_type]]$fac_0_longit_0$y$se[it] <- t_test_curr$stderr
      resl[[sim_type]]$fac_0_longit_0$y$pval[it] <- t_test_curr$p.value
    
      ##########################################################
      # Multiple time point analysis of single measurement
      treatment_means_form <- as.formula(paste0("y ~ 1 + arm:(", baseline_form_curr, ")"))
      contrasts_vs_placebo_form <- as.formula(paste0("y ~ -1 + time_fac + active_binary:(", baseline_form_curr, ")"))
      fixed_form_y <- switch(estimate_type,
                                treatment_means = treatment_means_form, 
                                contrasts_vs_placebo = contrasts_vs_placebo_form)
      lme_fac_0_longit_1 <- nlme::lme(fixed = fixed_form_y,
                                      random = ~ 1 | id_fac,
                                      data = results_df,
                                      correlation = nlme::corCAR1(form = ~ Time),
                                      control = nlme::lmeControl(returnObject = TRUE),
                                      method = "ML")
      summ_fac_0_longit_1 <- summary(lme_fac_0_longit_1)
      summ_fac_0_longit_1
      if (estimate_type == "contrasts_vs_placebo") {
        coef_name <- paste0("active_binary:time_", control$time_pt_eval, "_binary")
      }
      df_curr <- summ_fac_0_longit_1$tTable[coef_name, "DF"]
      resl[[sim_type]]$fac_0_longit_1$y$est[it] <- summ_fac_0_longit_1$tTable[coef_name, "Value"]
      resl[[sim_type]]$fac_0_longit_1$y$se[it] <- summ_fac_0_longit_1$tTable[coef_name, "Std.Error"]
      resl[[sim_type]]$fac_0_longit_1$y$pval[it] <- summ_fac_0_longit_1$tTable[coef_name, "p-value"]
      resl[[sim_type]]$fac_0_longit_1$y$low[it] <- resl[[sim_type]]$fac_0_longit_1$y$est[it] - 
        resl[[sim_type]]$fac_0_longit_1$y$se[it] * qt(p = 1 - control$ci_alpha / 2, df = df_curr)
      resl[[sim_type]]$fac_0_longit_1$y$upp[it] <- resl[[sim_type]]$fac_0_longit_1$y$est[it] + 
        resl[[sim_type]]$fac_0_longit_1$y$se[it] * qt(p = 1 - control$ci_alpha / 2, df = df_curr)
      fac_0_longit_1_pval <- summ_fac_0_longit_1$tTable[coef_name, "p-value"]
      
      ##########################################################
      # Single time point analysis of factors
      z_est_vec_at_estim_t<-z_se_vec_at_estim_t <- c()
      for (fac_curr in 1:n_fac) {
        z_nam <- paste0("z", fac_curr)
        active_z_single_t <- results_df[results_df$arm == "Active" & results_df$Time == control$time_pt_eval, z_nam]
        placebo_z_single_t <- results_df[results_df$arm == "Placebo" & results_df$Time == control$time_pt_eval, z_nam]
        t_test_curr <- t.test(x = active_z_single_t, y = placebo_z_single_t, var.equal = TRUE)
        df_curr <- t_test_curr$parameter
        resl[[sim_type]]$fac_1_longit_0[[z_nam]]$est[it] <- t_test_curr$estimate[1] - t_test_curr$estimate[2]
        resl[[sim_type]]$fac_1_longit_0[[z_nam]]$low[it] <- t_test_curr$conf.int[1]
        resl[[sim_type]]$fac_1_longit_0[[z_nam]]$upp[it] <- t_test_curr$conf.int[2]
        resl[[sim_type]]$fac_1_longit_0[[z_nam]]$se[it] <- t_test_curr$stderr
        resl[[sim_type]]$fac_1_longit_0[[z_nam]]$pval[it] <- t_test_curr$p.value
        z_est_vec_at_estim_t[fac_curr] <- resl[[sim_type]]$fac_1_longit_0[[z_nam]]$est[it]
        z_se_vec_at_estim_t[fac_curr] <- resl[[sim_type]]$fac_1_longit_0[[z_nam]]$se[it]
      }
      
      resl[[sim_type]]$fac_1_longit_0$y$est[it] <- t(mfd_object$A[y_nam, ]) %*% z_est_vec_at_estim_t * mfd_object$sd_scale_factor[y_nam]
      resl[[sim_type]]$fac_1_longit_0$y$se[it] <- sqrt(t(mfd_object$A[paste0("Y", control$y_ind_use), ]) %*% 
                                             diag((z_se_vec_at_estim_t)^2) %*% 
                                             mfd_object$A[paste0("Y", control$y_ind_use), ]) * mfd_object$sd_scale_factor[y_nam]
      resl[[sim_type]]$fac_1_longit_0$y$low[it] <- resl[[sim_type]]$fac_1_longit_0$y$est[it] - 
        resl[[sim_type]]$fac_1_longit_0$y$se[it] * qt(p = 1 - control$ci_alpha / 2, df = df_curr)
      resl[[sim_type]]$fac_1_longit_0$y$upp[it] <- resl[[sim_type]]$fac_1_longit_0$y$est[it] + 
        resl[[sim_type]]$fac_1_longit_0$y$se[it] * qt(p = 1 - control$ci_alpha / 2, df = df_curr)
      resl[[sim_type]]$fac_1_longit_0$y$pval[it] <- pt(q = abs(resl[[sim_type]]$fac_1_longit_0$y$est[it] / resl[[sim_type]]$fac_1_longit_0$y$se[it]), 
                                           df = df_curr, lower.tail = F) * 2
    
      ##########################################################
      # Multiple time point factor analysis
      z_treatment_means_form <- as.formula(paste0("z ~ 1 + arm:(", baseline_form_curr, ")"))
      z_contrasts_vs_placebo_form <- as.formula(paste0("z ~ -1 + time_fac + active_binary:(", baseline_form_curr, ")"))
      fixed_form_z <- switch(estimate_type,
                                  treatment_means = z_treatment_means_form,
                                  contrasts_vs_placebo = z_contrasts_vs_placebo_form)
      if (estimate_type == "contrasts_vs_placebo") {
        coef_name <- paste0("active_binary:time_", control$time_pt_eval, "_binary")
      }
      z_est_vec_at_estim_t<-z_se_vec_at_estim_t <- c()
      for (fac_curr in 1:n_fac) {
        z_nam <- paste0("z", fac_curr)
        s_nam <- paste0("s", fac_curr)
        results_df$z <- results_df[, z_nam]
        results_df$s <- results_df[, s_nam]
        
        lme_fac_1_longit_1 <- nlme::lme(fixed = fixed_form_z,
                                      random = ~ 1 | id_fac,
                                      data = results_df,
                                      correlation = nlme::corCAR1(form = ~ Time),
                                      control = nlme::lmeControl(returnObject = TRUE),
                                      method = "ML")
        summ_fac_1_longit_1 <- summary(lme_fac_1_longit_1)
        summ_fac_1_longit_1
        df_curr <- summ_fac_1_longit_1$tTable[coef_name, "DF"]
        resl[[sim_type]]$fac_1_longit_1[[z_nam]]$est[it] <- summ_fac_1_longit_1$tTable[coef_name, "Value"]
        resl[[sim_type]]$fac_1_longit_1[[z_nam]]$se[it] <- summ_fac_1_longit_1$tTable[coef_name, "Std.Error"]
        resl[[sim_type]]$fac_1_longit_1[[z_nam]]$pval[it] <- summ_fac_1_longit_1$tTable[coef_name, "p-value"]
        fac_1_longit_1_pval <- summ_fac_1_longit_1$tTable[coef_name, "p-value"]
        resl[[sim_type]]$fac_1_longit_1[[z_nam]]$low[it] <- resl[[sim_type]]$fac_1_longit_1[[z_nam]]$est[it] - 
          resl[[sim_type]]$fac_1_longit_1[[z_nam]]$se[it] * qt(p = 1 - control$ci_alpha / 2, df = df_curr)
        resl[[sim_type]]$fac_1_longit_1[[z_nam]]$upp[it] <- resl[[sim_type]]$fac_1_longit_1[[z_nam]]$est[it] +
          resl[[sim_type]]$fac_1_longit_1[[z_nam]]$se[it] * qt(p = 1 - control$ci_alpha / 2, df = df_curr)
        z_est_vec_at_estim_t[fac_curr] <- resl[[sim_type]]$fac_1_longit_1[[z_nam]]$est[it]
        z_se_vec_at_estim_t[fac_curr] <- resl[[sim_type]]$fac_1_longit_1[[z_nam]]$se[it]
      }
      resl[[sim_type]]$fac_1_longit_1$y$est[it] <- t(mfd_object$A[y_nam, ]) %*% z_est_vec_at_estim_t * mfd_object$sd_scale_factor[y_nam]
      resl[[sim_type]]$fac_1_longit_1$y$se[it] <- sqrt(sum(mfd_object$A[paste0("Y", control$y_ind_use), ]^2 * z_se_vec_at_estim_t^2)) * 
        mfd_object$sd_scale_factor[y_nam]
      resl[[sim_type]]$fac_1_longit_1$y$low[it] <- resl[[sim_type]]$fac_1_longit_1$y$est[it] - 
        resl[[sim_type]]$fac_1_longit_1$y$se[it] * qt(p = 1 - control$ci_alpha / 2, df = df_curr)
      resl[[sim_type]]$fac_1_longit_1$y$upp[it] <- resl[[sim_type]]$fac_1_longit_1$y$est[it] + 
        resl[[sim_type]]$fac_1_longit_1$y$se[it] * qt(p = 1 - control$ci_alpha / 2, df = df_curr)
      resl[[sim_type]]$fac_1_longit_1$y$pval[it] <- pt(q = abs(resl[[sim_type]]$fac_1_longit_1$y$est[it] / resl[[sim_type]]$fac_1_longit_1$y$se[it]), 
                                           df = df_curr, lower.tail = F) * 2
      
    }
  }      
  
  resl[["non_null"]]$true_estimand <- sim_non_null$E_Y[as.character(control$time_pt_eval), y_nam, "trt_1"] - 
    sim_non_null$E_Y[as.character(control$time_pt_eval), y_nam, "trt_0"]
  resl[["null"]]$true_estimand <- 0
  resl[["non_null"]]$example_sim <- sim_non_null
  resl[["null"]]$example_sim <- sim_null
  saveRDS(resl, file = file.path("output", 
                                 paste0("resl_parameter_set_", TASK_ID, ".RDS")))
}
