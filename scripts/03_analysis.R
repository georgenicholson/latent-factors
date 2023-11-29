
set_baseline_arm_effects_equal <- TRUE
control <- get_control_parameters_analysis()

analysis_names <- names(control$do_list)

##########################################################################################
# Gather study data subsets
analysis_sub_data <- list()
for (analysis_curr in analysis_names) {
  d_time_sub <- d_time %>% 
    filter(study %in% control$do_list[[analysis_curr]]$stage1) %>%
    filter(USUBJID %in% d_subject$USUBJID) %>%
    filter(AVISITN %in% c(0, 2, 4, 8, 12, 16))
  d_subject_sub <- d_subject %>% 
    filter(USUBJID %in% d_time_sub$USUBJID)
  ###################################################
  # Normal score transformation happening here
  d_time_sub <- mutate_if(as.data.frame(d_time_sub), colnames(d_time_sub) %in% nms_params_v2$PARAMCD, normal_quantile_map)
  colnames(d_time_sub)[colnames(d_time_sub) %in% nms_params_v2$PARAMCD] <- 
    nms_params_v2$PARAM[match(colnames(d_time_sub)[colnames(d_time_sub) %in% nms_params_v2$PARAMCD], nms_params_v2$PARAMCD)]
  d_time_sub <- as_tibble(d_time_sub)  
  analysis_sub_data[[analysis_curr]] <- list(time = d_time_sub, subject = d_subject_sub)
  analysis_sub_data[[analysis_curr]]$subject$newTRT <- gsub(" ", "", analysis_sub_data[[analysis_curr]]$subject$newTRT)
  
}

examine_scale_of_qn_outputs <- FALSE
if (examine_scale_of_qn_outputs) {
  par(mfrow = c(3, 4))
  for (parc in intersect(nms_params_v2$PARAM, names(analysis_sub_data[[analysis_curr]]$time))) {
    yc <- na.omit(as.data.frame(analysis_sub_data[[analysis_curr]]$time)[, parc])
    hist(yc)
    print(t(c(mean(yc), sd(yc))))
  }
}

##########################################################################################
# Stage 1 fits
mfd_object_list <- list()
for (analysis_curr in analysis_names) {
  mfd_object_list[[analysis_curr]] <- mlfa::mfd(data = analysis_sub_data[[analysis_curr]]$time, 
                                                vars = nms_params_v2$PARAM, 
                                                loadings.n_PCs = control$n_fac, 
                                                seed = 1,
                                                time = "AVISITN", 
                                                subject = "USUBJID",
                                                normalize.type = c("per_meas_standard", "standardize", "per_center_baseline_scale")[1],
                                                Time_base = 0, 
                                                standardize_A_type = c("none", "L2", "largest_abs")[3],
                                                loadings.sparse_abs_thres = 0,
                                                loadings.rotation = c("varimax", "promax")[1],
                                                scores.type = 'a_posteriori')
  
}


##########################################################################################
# Ensure the loadings matrices columns correspond across studies
analysis_to_match_loadings_to <- 'PsA_F2312'
for (analysis_curr in setdiff(analysis_names, analysis_to_match_loadings_to)) {
  prox_mat <- t(mfd_object_list[[analysis_curr]]$A) %*% mfd_object_list[[analysis_to_match_loadings_to]]$A
  index_move_from <- apply(prox_mat, 2, which.max)
  mfd_object_list[[analysis_curr]]$A <- mfd_object_list[[analysis_curr]]$A[, index_move_from]
  mfd_object_list[[analysis_curr]]$Z <- mfd_object_list[[analysis_curr]]$Z[, index_move_from]
  mfd_object_list[[analysis_curr]]$S <- mfd_object_list[[analysis_curr]]$S[, index_move_from]
  mfd_object_list[[analysis_curr]]$data_out[, paste0("z", 1:3)] <- 
    mfd_object_list[[analysis_curr]]$data_out[, paste0("z", 1:3)[index_move_from]]
  mfd_object_list[[analysis_curr]]$data_out[, paste0("s", 1:3)] <- 
    mfd_object_list[[analysis_curr]]$data_out[, paste0("s", 1:3)[index_move_from]]
}

########################################
# Use lme() to fit trajectories
analysis_results_list <- list()
vars_for_sim <- data.frame()
analysis_curr <- analysis_names[2]
for (analysis_curr in analysis_names) {
  placebo_label <- "Placebo"
  multiple_studies_flag <- length(control$do_list[[analysis_curr]]$stage2) > 1
  studies_in_stage_2 <- control$do_list[[analysis_curr]]$stage2
  results_df <- as.data.frame(mfd_object_list[[analysis_curr]]$data_out)
  results_df$STUDYID <- sapply(strsplit(results_df$USUBJID, split = "_"), function(x) x[[1]][1])
  results_df$study <- gsub("CAIN457", "", results_df$STUDYID)
  results_df <- results_df[results_df$study %in% studies_in_stage_2, ]
  if (multiple_studies_flag) {
    results_df$study <- C(factor(results_df$study), contr = contr.sum) # sum-to-zero contrasts for study
  }  
  results_df$sex <- d_subject[match(results_df$USUBJID, d_subject$USUBJID), ] %>% pull(SEX)
  results_df$sex <- C(factor(results_df$sex), contr = contr.sum) # sum-to-zero contrasts for sex
  results_df$arm <- as.factor(analysis_sub_data[[analysis_curr]]$subject[
    match(results_df$USUBJID, analysis_sub_data[[analysis_curr]]$subject$USUBJID), ] %>% pull(newTRT))
  results_df$arm <- relevel(results_df$arm, ref = placebo_label)
  results_df$time_fac <- relevel(factor(results_df$AVISITN), ref = "0")
  results_df$resid_grp <- as.factor(1:nrow(results_df))
  t_unique <- sort(unique(results_df$AVISITN))
  nt <- length(t_unique)
  arm_all <- unique(results_df$arm)
  n_arm <- length(arm_all)
  arm_no_plac <- setdiff(arm_all, placebo_label)
  for (time_curr in t_unique) {
    results_df[, paste0("time_", time_curr, "_binary")] <- ifelse(results_df$AVISITN == time_curr, 1, 0)
  }
  for (arm_curr in arm_all) {
    arm_curr_bin_name <- paste0("arm_", arm_curr, "_binary")
    results_df[, arm_curr_bin_name] <- ifelse(results_df$arm == arm_curr, 1, 0)
  }    
  m_coef_mat<-u_coef_mat<-l_coef_mat <- matrix(NA, nt, n_arm, dimnames = list(t_unique, arm_all))
  results_list <- list()
  fac_curr <- 1
  for (fac_curr in 1:control$n_fac) {
    results_list[[fac_curr]] <- list()
    results_df$z <- results_df[, paste0("z", fac_curr)]
    results_df$s <- results_df[, paste0("s", fac_curr)]
    estimate_type <- c("treatment_means", "contrasts_vs_placebo")[2]
    for (estimate_type in control$estimate_types) {
      model_formula_object <- get_formula(estimate_type = estimate_type, 
                                          set_baseline_arm_effects_equal = set_baseline_arm_effects_equal, 
                                          multiple_studies_flag = multiple_studies_flag, 
                                          arm_all = arm_all, 
                                          t_unique = t_unique, 
                                          placebo_label = placebo_label)
      fixed_form_curr <- model_formula_object$formula
      
      coef_names <- model_formula_object$coef_names
      lme_out <- nlme::lme(fixed = fixed_form_curr,
                           random = ~ 1 | USUBJID, 
                           data = results_df,
                           correlation = nlme::corCAR1(form = ~ AVISITN),
                           control = nlme::lmeControl(returnObject = TRUE),
                           method = "ML")
      summ_curr <- summary(lme_out)
      var_tab <- nlme::VarCorr(summ_curr)
      cis <- nlme::intervals(lme_out)
      if (estimate_type == "contrasts_vs_placebo") {
        pick_index <- coef_names[["150mg"]][grep("16", coef_names[["150mg"]])]
        add <- data.frame(analysis = analysis_curr,
                          fac = fac_curr,
                          estimate_type = estimate_type,
                          RI_sd_est = cis$reStruct$USUBJID$est., 
                          RI_sd_low = cis$reStruct$USUBJID$lower, 
                          RI_sd_upp = cis$reStruct$USUBJID$upper, 
                          z_resid_sd_est = cis$sigma["est."], 
                          z_resid_sd_low = cis$sigma["lower"], 
                          z_resid_sd_upp = cis$sigma["upper"], 
                          Phi_est = nlme::intervals(lme_out)$corStruct["Phi", "est."],
                          Phi_low = nlme::intervals(lme_out)$corStruct["Phi", "lower"],
                          Phi_est = nlme::intervals(lme_out)$corStruct["Phi", "upper"],
                          week_16_150mg_est = cis$fixed[pick_index, "est."],
                          week_16_150mg_low = cis$fixed[pick_index, "lower"],
                          week_16_150mg_upp = cis$fixed[pick_index, "upper"],
                          y_resid_sd_est = sqrt(mfd_object_list[[analysis_curr]]$sigma_square_ml))
        vars_for_sim <- rbind(vars_for_sim, add)
      }
      cis
      mn_fixed <- summ_curr$coefficients$fixed
      cov_fixed <- summ_curr$varFix
      se_fixed <- sqrt(diag(cov_fixed))
      cis$fixed <- rbind(cis$fixed, matrix(0, ncol = 3, nrow = 1, dimnames = list("zero", colnames(cis$fixed))))
      for (arm_curr in arm_all) {
        m_coef_mat[, arm_curr] <- cis$fixed[coef_names[[arm_curr]], "est."]
        u_coef_mat[, arm_curr] <- cis$fixed[coef_names[[arm_curr]], "upper"]
        l_coef_mat[, arm_curr] <- cis$fixed[coef_names[[arm_curr]], "lower"]
      }
      results_list[[fac_curr]][[estimate_type]] <- list(m = m_coef_mat,
                                                        u = u_coef_mat,
                                                        l = l_coef_mat)
    }
  }
  analysis_results_list[[analysis_curr]] <- results_list
  analysis_results_list[[analysis_curr]]$mfd_object <- mfd_object_list[[analysis_curr]]
  analysis_results_list[[analysis_curr]]$arm_all <- arm_all
  analysis_results_list[[analysis_curr]]$arm_no_plac <- arm_no_plac
  
}

vars_for_sim$RI_var <- as.numeric(vars_for_sim$RI_sd_est^2)
vars_for_sim$z_resid_var <- as.numeric(vars_for_sim$z_resid_sd_est^2)
vars_for_sim$y_resid_var <- as.numeric(vars_for_sim$y_resid_sd_est^2)
vars_for_sim$var_tot <- vars_for_sim$RI_var + vars_for_sim$z_resid_var + vars_for_sim$y_resid_var

saveRDS(analysis_results_list, file = file.path("output", "analysis_results_list.RDS"))
saveRDS(vars_for_sim, file = file.path("output", "vars_for_sim.RDS"))



