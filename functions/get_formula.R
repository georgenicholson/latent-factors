get_formula <- function(estimate_type, 
                        set_baseline_arm_effects_equal, 
                        multiple_studies_flag, 
                        arm_all, 
                        t_unique, 
                        placebo_label = "Placebo") {
  coef_names <- list()
  if (!"0" %in% t_unique){
    stop("Baseline week 0 not present in t_unique")
  }
  t_unique_no_baseline <- setdiff(t_unique, "0")
  sex_study_string <- ifelse(multiple_studies_flag, "sex + study", "sex")
  add_times_wo_baseline <- paste0("(", paste(paste0("time_", setdiff(t_unique, 0), "_binary"), collapse = " + "), ")")
  add_times_with_baseline <- paste0("(", paste(paste0("time_", t_unique, "_binary"), collapse = " + "), ")")
  add_arms_no_placebo <- paste0("(", paste(paste0("arm_", setdiff(arm_all, placebo_label), "_binary"), collapse = " + "), ")")
  if (estimate_type == "treatment_means") {
    if (set_baseline_arm_effects_equal) {
      form_out_char <- paste0("z ~ 1 + ", sex_study_string, " + ", add_times_wo_baseline, ":(arm + ", sex_study_string, ")")
      for (arm_curr in arm_all) {
        coef_names[[arm_curr]] <- c("(Intercept)", paste0("time_", t_unique_no_baseline, "_binary:arm", arm_curr))
      }
    } else {
      form_out_char <- paste0("z ~ -1 + time_fac:arm + time_fac:(", sex_study_string, ")")
      for (arm_curr in arm_all) {
        coef_names[[arm_curr]] <- c(paste0("time_fac", t_unique, ":arm", arm_curr))
      }
    }
  }
  if (estimate_type == "contrasts_vs_placebo") {
    if (set_baseline_arm_effects_equal) {
      form_out_char <- paste0("z ~ time_fac + ", 
                              add_times_wo_baseline, ":", add_arms_no_placebo, " + ", 
                              "time_fac:(", sex_study_string, ")")
      fixed_form_curr <- as.formula(form_out_char)
      X <- model.matrix(fixed_form_curr, model.frame(fixed_form_curr, data = results_df))
      qr(X)$rank
      dim(X)
      colnames(X)
      for (arm_curr in arm_all) {
        if (arm_curr == placebo_label) {
          coef_names[[arm_curr]] <- c("(Intercept)", paste0("time_fac", t_unique_no_baseline))
        } else {
          coef_names[[arm_curr]] <- c("zero", paste0("time_", t_unique_no_baseline, "_binary:arm_", arm_curr, "_binary"))
        }
      }
      
    } else {
      form_out_char <- paste0("z ~ -1 + time_fac + time_fac:arm + time_fac:(", sex_study_string, ")")
      for (arm_curr in arm_all) {
        if (arm_curr == placebo_label) {
          coef_names[[arm_curr]] <- paste0("time_fac", t_unique)
        } else {
          coef_names[[arm_curr]] <- paste0("time_fac", t_unique, ":arm", arm_curr)
        }
      }
    }
  }
  form_out <- as.formula(form_out_char)
  return(list(formula = form_out, coef_names = coef_names))
}

