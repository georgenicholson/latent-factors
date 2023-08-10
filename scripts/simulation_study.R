# ##################################################################
# # Process command line arguments
# ##################################################################
# source("functions/process_command_line_args.R")
# args <- process_command_line_args(args_default = list(TASK_ID = 500,
#                                                       N_TASKS = 1000,
#                                                       CONTROL_NUMBER = 1,
#                                                       SNP_NUM = 8))
# TASK_ID <- args$TASK_ID
# N_TASKS <- args$N_TASKS
# CONTROL_NUMBER <- args$CONTROL_NUMBER
# SNP_NUM <- args$SNP_NUM
# snp_nam_ord <- readRDS(file = "output/linear_model_fitting/snp_nam_ord.RDS")
# snp_fit <- snp_nam_ord[SNP_NUM]

# 
# devtools::load_all("mlfa_package")  # not loading it, because we are currently working with a development script for mfd which is loaded in 03_analysis.R, see source("scripts/mfd_develop.R")
# devtools::install("mlfa_package")
# library(mlfa)


#######################################
# Source functions
fn_files_to_source <- list.files("functions/", full.names = TRUE)
for (file_curr in fn_files_to_source) {
  source(file_curr)
}

# #######################################
# # Folder structure
# DIR <- get_DIR()

control <- get_control_parameters()
analysis_names <- names(control$do_list)
placebo_label <- "Placebo"
estimate_type <- c("treatment_means", "contrasts_vs_placebo")[2]


# Varying params

#set param's by default at something resembling data
# then vary each set of param's away from there individually
n_subj <- c(50)
# These parameters determine the autocorrelated residual component of z
theta_resid_sd <- 1
cor_ar1_lag <- c(1 - 1e-5, .99, .975, .95, .75, .5, 0)[1]
# This parameter determine the random intercept component of z
theta_RI_sd <- 1
# This parameter determine the IID residual noise in y
y_resid_sd <- c(0.01, .25, 1, 2, 10)[4]
endpoint_proportions <- list(c(0.6, 0.3, 0.1),
                             c(1, 1, 1, 1) / 4)[[1]]
missing_proportion <- c(0, .25)[1]
beta_Emax <- c(0, .1, .25, .5)[3]
set_baseline_arm_effects_equal <- FALSE

# Fixed params
TRT_prob <- 0.5
E0_null <- 0
E50_null <- 7
Emax_null <- 3
beta_E0 <- 0
beta_E50 <- 0
n_fac <- length(endpoint_proportions)
n_var <- 12
weeks_measure <- c(0, 1, 2, 4, 8, 12, 16)
time_pt_eval <- "16"
var_names <- paste0("Y", 1:n_var)
n_its <- 100
ci_alpha <- 0.05
pval_thresh <- 0.05
y_ind_use <- 1

model_names <- c("fac_0_longit_0", "fac_0_longit_1", "fac_1_longit_0", "fac_1_longit_1")
restab <- data.frame(it = 1:n_its, est = NA, low = NA, upp = NA, se = NA, pval = NA)
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
pval_df <- data.frame()
for (it in 1:n_its) {
# for (it in 1) {
    if (it %% 25 == 0) {
    print(it)
  }
  for (sim_type in c("null", "non_null")) {
    ##########################################
    # Simulate data 
    sim_curr <- simulate_data(seed = it,
                    n = n_subj,
                    p = n_var,
                    k = n_fac,
                    TRT_prob = TRT_prob,
                    endpoint_proportions = endpoint_proportions,
                    time = weeks_measure,
                    cor_ar1_lag = cor_ar1_lag,
                    E0_null = E0_null,
                    E50_null = E50_null,
                    Emax_null = Emax_null,
                    beta_E0 = beta_E0,
                    beta_E50 = beta_E50,
                    beta_Emax = switch(sim_type,
                                       non_null = beta_Emax,
                                       null = 0),
                    theta_RI_sd = theta_RI_sd,
                    theta_resid_sd = theta_resid_sd,
                    y_resid_sd = y_resid_sd)
  
    if (sim_type == "non_null") {
      sim_non_null <- sim_curr
    }
    data_in <- as_tibble(sim_curr$Y)
    sub_data <- list(time = as_tibble(sim_curr$Y),
                     subject = as_tibble(sim_curr$dbase))
    
    
    y_test <- as.matrix(sim_curr$Z[, paste0("Z", 1:3)]) %*% sim_curr$A_loadings
    # par(mfrow = c(1, 2))
    # plot(y_test[, 1], sim_curr$Y[, "Y1"])
    # abline(0, 1)
    # plot(y_test[, 1], sim_curr$Y[, "Y1"] + rnorm(n = nrow(y_test), sd = 2))
    # abline(0, 1)
    ##########################################
    # Stage 1: estimate sparse loadings A and scores z
    mfd_object <- mlfa::mfd(data = sub_data$time, 
                vars = var_names, 
                loadings.n_PCs = n_fac, 
                seed = 1,
                time = "Time", 
                subject = "ID",
                normalize.type = c("per_meas_standard", "standardize", "per_center_baseline_scale")[1],
                Time_base = 0, 
                standardize_A_type = c("none", "L2", "largest_abs")[3],
                loadings.sparse_abs_thres = 0,
                loadings.rotation = c("varimax", "promax")[1],
                scores.type = 'a_posteriori', 
                trace = FALSE)
    # par(mfrow = c(1, 2))
    # image(t(mfd_object$A))
    # image(sim_curr$A_loadings)
  
    # apply(mfd_object$Z, 2, sd)
    # apply(mfd_object$Z, 2, mean)
    # apply(sim_curr$Z[, paste0("Z", 1:3)], 2, sd)
    # apply(as.matrix(sim_curr$Z[, paste0("Z", 1:3)]) %*% as.matrix(sim_curr$A_loadings), 2, sd)
    # apply(as.matrix(sim_curr$Y[, paste0("Y", 1:n_var)]), 2, sd)
    # apply(mfd_object$Z %*% t(as.matrix(mfd_object$A)) * mfd_object$sd_scale_factor, 2, sd)
     
    ##########################################
    # Assemble data for Stage 2
    arm_all <- unique(sub_data$subject$TRT)
    n_arm <- length(arm_all)
    results_df <- as.data.frame(mfd_object$data_out)
    y_nam <- paste0("Y", y_ind_use)
    z_ind_use <- which.max(c(cor(sub_data$time %>% pull(paste0("Y", y_ind_use)), results_df[, paste0("z", 1:n_fac)])))
    results_df$arm <- sub_data$subject$TRT[match(results_df$ID, sub_data$subject$ID)]
    results_df$arm <- relevel(factor(results_df$arm), ref = "Placebo")
    results_df$time_fac <- factor(results_df$Time)
    results_df$id_fac <- factor(results_df$ID)
    nt <- length(weeks_measure)
    m_coef_mat<-u_coef_mat<-l_coef_mat <- matrix(NA, nt, n_arm, dimnames = list(weeks_measure, arm_all))
    results_list <- list()
    results_df$y <- sub_data$time %>% pull(paste0("Y", y_ind_use))
    # fac_curr <- 1
    # results_list[[fac_curr]] <- list()
    # results_df$z <- results_df[, paste0("z", fac_curr)]
    # results_df$s <- results_df[, paste0("s", fac_curr)]
    
    
    
    # #################
    # plot(as.matrix(results_df[, c("z1", "z2", "z3")]) %*% sim_curr$A_loadings[, paste0("P", y_ind_use)] ,
    #      results_df$y * mfd_object$sd_scale_factor[y_ind_use])
    # 
    # plot(as.matrix(results_df[, c("z1", "z2", "z3")]) %*% mfd_object$A[paste0("Y", y_ind_use), ] ,
    #      results_df$y * mfd_object$sd_scale_factor[y_ind_use])
    # abline(0, 1)
    # sd(results_df$y)
    # apply(as.matrix(results_df[, c("z1", "z2", "z3")]), 2, sd)
    # apply(sim_data$time, 2, sd)
    # ##################
    for (time_curr in weeks_measure) {
      results_df[, paste0("time_", time_curr, "_binary")] <- ifelse(results_df$Time == time_curr, 1, 0)
    }
    non_baseline_form <- paste(paste0("time_", setdiff(weeks_measure, 0), "_binary"), collapse = " + ")
    include_baseline_form <- paste(paste0("time_", weeks_measure, "_binary"), collapse = " + ")
    if (set_baseline_arm_effects_equal) {
      baseline_form_curr <- non_baseline_form
    } else {
      baseline_form_curr <- include_baseline_form
    }
    results_df$active_binary <- ifelse(results_df$arm == "Active", 1, 0)
    
    ##########################################################
    # Single time point analysis of single measurement
    active_y_single_t <- results_df[results_df$arm == "Active" & results_df$Time == time_pt_eval, "y"]
    placebo_y_single_t <- results_df[results_df$arm == "Placebo" & results_df$Time == time_pt_eval, "y"]
    t_test_curr <- t.test(x = active_y_single_t, y = placebo_y_single_t)
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
    
    # lme_fac_0_longit_1 <- nlme::lme(fixed = fixed_form_y,
    #                      random = ~ 1 | id_fac,
    #                      data = results_df,
    #                      control = nlme::lmeControl(returnObject = TRUE),
    #                      method = "ML")
    lme_fac_0_longit_1 <- nlme::lme(fixed = fixed_form_y,
                                    random = ~ 1 | id_fac,
                                    data = results_df,
                                    correlation = nlme::corCAR1(form = ~ Time),
                                    control = nlme::lmeControl(returnObject = TRUE),
                                    method = "ML")
    summ_fac_0_longit_1 <- summary(lme_fac_0_longit_1)
    summ_fac_0_longit_1
    if (estimate_type == "contrasts_vs_placebo") {
      coef_name <- paste0("active_binary:time_", time_pt_eval, "_binary")
    }
    df_curr <- summ_fac_0_longit_1$tTable[coef_name, "DF"]
    resl[[sim_type]]$fac_0_longit_1$y$est[it] <- summ_fac_0_longit_1$tTable[coef_name, "Value"]
    resl[[sim_type]]$fac_0_longit_1$y$se[it] <- summ_fac_0_longit_1$tTable[coef_name, "Std.Error"]
    resl[[sim_type]]$fac_0_longit_1$y$pval[it] <- summ_fac_0_longit_1$tTable[coef_name, "p-value"]
    resl[[sim_type]]$fac_0_longit_1$y$low[it] <- resl[[sim_type]]$fac_0_longit_1$y$est[it] - 
      resl[[sim_type]]$fac_0_longit_1$y$se[it] * qt(p = 1 - ci_alpha / 2, df = df_curr)
    resl[[sim_type]]$fac_0_longit_1$y$upp[it] <- resl[[sim_type]]$fac_0_longit_1$y$est[it] + 
      resl[[sim_type]]$fac_0_longit_1$y$se[it] * qt(p = 1 - ci_alpha / 2, df = df_curr)
    
    fac_0_longit_1_pval <- summ_fac_0_longit_1$tTable[coef_name, "p-value"]
  
    
    ##########################################################
    # Single time point analysis of factors
    z_est_vec_at_estim_t<-z_se_vec_at_estim_t <- c()
    for (fac_curr in 1:n_fac) {
      z_nam <- paste0("z", fac_curr)
      active_z_single_t <- results_df[results_df$arm == "Active" & results_df$Time == time_pt_eval, z_nam]
      placebo_z_single_t <- results_df[results_df$arm == "Placebo" & results_df$Time == time_pt_eval, z_nam]
      t_test_curr <- t.test(x = active_z_single_t, y = placebo_z_single_t)
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
    resl[[sim_type]]$fac_1_longit_0$y$se[it] <- sqrt(t(mfd_object$A[paste0("Y", y_ind_use), ]) %*% 
                                           diag((z_se_vec_at_estim_t)^2) %*% 
                                           mfd_object$A[paste0("Y", y_ind_use), ]) * mfd_object$sd_scale_factor[y_nam]
    resl[[sim_type]]$fac_1_longit_0$y$low[it] <- resl[[sim_type]]$fac_1_longit_0$y$est[it] - 
      resl[[sim_type]]$fac_1_longit_0$y$se[it] * qt(p = 1 - ci_alpha / 2, df = df_curr)
    resl[[sim_type]]$fac_1_longit_0$y$upp[it] <- resl[[sim_type]]$fac_1_longit_0$y$est[it] + 
      resl[[sim_type]]$fac_1_longit_0$y$se[it] * qt(p = 1 - ci_alpha / 2, df = df_curr)
    resl[[sim_type]]$fac_1_longit_0$y$pval[it] <- pt(q = resl[[sim_type]]$fac_1_longit_0$y$est[it] / resl[[sim_type]]$fac_1_longit_0$y$se[it], 
                                         df = df_curr, lower.tail = F) * 2
  
    ##########################################################
    # Multiple time point factor analysis
    # z_treatment_means_form <- as.formula(paste0("z ~ 1 + arm:(", baseline_form_curr, ")"))
    # z_contrasts_vs_placebo_form <- as.formula(paste0("z ~ -1 + time_fac + active_binary:(", baseline_form_curr, ")"))
    z_treatment_means_form <- as.formula(paste0("z ~ 1 + arm:(", baseline_form_curr, ")"))
    z_contrasts_vs_placebo_form <- as.formula(paste0("z ~ -1 + time_fac + active_binary:(", baseline_form_curr, ")"))
    # z_contrasts_vs_placebo_form <- as.formula(paste0("z ~ -1 + time_fac + arm:time_fac"))
    fixed_form_z <- switch(estimate_type,
                                treatment_means = z_treatment_means_form,
                                contrasts_vs_placebo = z_contrasts_vs_placebo_form)
    if (estimate_type == "contrasts_vs_placebo") {
      coef_name <- paste0("active_binary:time_", time_pt_eval, "_binary")
    }
    z_est_vec_at_estim_t<-z_se_vec_at_estim_t <- c()
    for (fac_curr in 1:n_fac) {
      z_nam <- paste0("z", fac_curr)
      s_nam <- paste0("s", fac_curr)
      results_df$z <- results_df[, z_nam]
      results_df$s <- results_df[, s_nam]
      
      
      # lm(fixed_form_z, data = results_df)
      # mean(results_df$z[results_df$Time == 16]) - mean(results_df$z[results_df$Time == 0])
      lme_fac_1_longit_1 <- nlme::lme(fixed = fixed_form_z,
                                    random = ~ 1 | id_fac,
                                    # weights = nlme::varFixed(~ s^2),
                                    data = results_df,
                                    # correlation = nlme::corAR1(form = ~ Time),
                                    correlation = nlme::corCAR1(form = ~ Time),
                                    control = nlme::lmeControl(returnObject = TRUE),
                                    # control = nlme::lmeControl(sigma = 1, returnObject = TRUE),
                                    method = "ML")
      summ_fac_1_longit_1 <- summary(lme_fac_1_longit_1)
      summ_fac_1_longit_1
      df_curr <- summ_fac_1_longit_1$tTable[coef_name, "DF"]
      resl[[sim_type]]$fac_1_longit_1[[z_nam]]$est[it] <- summ_fac_1_longit_1$tTable[coef_name, "Value"]
      resl[[sim_type]]$fac_1_longit_1[[z_nam]]$se[it] <- summ_fac_1_longit_1$tTable[coef_name, "Std.Error"]
      resl[[sim_type]]$fac_1_longit_1[[z_nam]]$pval[it] <- summ_fac_1_longit_1$tTable[coef_name, "p-value"]
      fac_1_longit_1_pval <- summ_fac_1_longit_1$tTable[coef_name, "p-value"]
      resl[[sim_type]]$fac_1_longit_1[[z_nam]]$low[it] <- resl[[sim_type]]$fac_1_longit_1[[z_nam]]$est[it] - 
        resl[[sim_type]]$fac_1_longit_1[[z_nam]]$se[it] * qt(p = 1 - ci_alpha / 2, df = df_curr)
      resl[[sim_type]]$fac_1_longit_1[[z_nam]]$upp[it] <- resl[[sim_type]]$fac_1_longit_1[[z_nam]]$est[it] +
        resl[[sim_type]]$fac_1_longit_1[[z_nam]]$se[it] * qt(p = 1 - ci_alpha / 2, df = df_curr)
      z_est_vec_at_estim_t[fac_curr] <- resl[[sim_type]]$fac_1_longit_1[[z_nam]]$est[it]
      z_se_vec_at_estim_t[fac_curr] <- resl[[sim_type]]$fac_1_longit_1[[z_nam]]$se[it]
    }
    resl[[sim_type]]$fac_1_longit_1$y$est[it] <- t(mfd_object$A[y_nam, ]) %*% z_est_vec_at_estim_t * mfd_object$sd_scale_factor[y_nam]
    resl[[sim_type]]$fac_1_longit_1$y$se[it] <- sqrt(sum(mfd_object$A[paste0("Y", y_ind_use), ]^2 * z_se_vec_at_estim_t^2)) * 
      mfd_object$sd_scale_factor[y_nam]
    resl[[sim_type]]$fac_1_longit_1$y$low[it] <- resl[[sim_type]]$fac_1_longit_1$y$est[it] - 
      resl[[sim_type]]$fac_1_longit_1$y$se[it] * qt(p = 1 - ci_alpha / 2, df = df_curr)
    resl[[sim_type]]$fac_1_longit_1$y$upp[it] <- resl[[sim_type]]$fac_1_longit_1$y$est[it] + 
      resl[[sim_type]]$fac_1_longit_1$y$se[it] * qt(p = 1 - ci_alpha / 2, df = df_curr)
    resl[[sim_type]]$fac_1_longit_1$y$pval[it] <- pt(q = resl[[sim_type]]$fac_1_longit_1$y$est[it] / resl[[sim_type]]$fac_1_longit_1$y$se[it], 
                                         df = df_curr, lower.tail = F) * 2
    
    # par(mfrow = c(1, 2))
    # y_df <- results_df
    # ymat_pl <- NULL
    # for (idc in 1:10) {
    #   ymat_pl <- cbind(ymat_pl, y_df[y_df$ID == idc, "y"])
    # }
    # matplot(ymat_pl, type = "l")
    # z_df <- results_df
    # zmat_pl <- NULL
    # for (idc in 1:10) {
    #   zmat_pl <- cbind(zmat_pl, z_df[z_df$ID == idc, "z"])
    # }
    # matplot(zmat_pl, type = "l")
    
    
  }
}      

# print(lapply(resl[["non_null"]], function(x) x$y[1:10, ]))


true_val <- sim_non_null$E_Y[time_pt_eval, y_nam, "trt_1"] - sim_curr$E_Y[time_pt_eval, y_nam, "trt_0"]
mae <- sapply(resl[["non_null"]], function(x) mean(abs(x$y$est - true_val)))
mse <- sapply(resl[["non_null"]], function(x) mean((x$y$est - true_val)^2))
bias_est <- sapply(resl[["non_null"]], function(x) mean(x$y$est - true_val))
bias_se <- sapply(resl[["non_null"]], function(x) sd(x$y$est - true_val) / sqrt(n_its))
sd_est <- sapply(resl[["non_null"]], function(x) sd(x$y$est - true_val))
# sd_se <- sapply(resl[["non_null"]], function(x) sd(x$y$est - true_val))
power <- sapply(resl[["non_null"]], function(x) mean(x$y$pval < pval_thresh))
power_se <- sapply(resl[["non_null"]], function(x) sd(x$y$pval < pval_thresh) / sqrt(n_its))
type1 <- sapply(resl[["null"]], function(x) mean((x$y$pval < pval_thresh)))
coverage <- sapply(resl[["non_null"]], function(x) mean(x$y$low <= true_val & x$y$upp >= true_val))
coverage_se <- sapply(resl[["non_null"]], function(x) sd(x$y$low <= true_val & x$y$upp >= true_val) / sqrt(n_its))
mn_ci_width <- sapply(resl[["non_null"]], function(x) mean(x$y$upp - x$y$low))

ndp <- 2
ndp_bias <- 3
ndp_mse <- 3
ndp_se <- 3
tab_res <- data.frame(meth = names(mse),
                      mae = round(mae, ndp),
                      mse = round(mse, ndp_mse),
                      bias = round(bias_est, ndp_bias),
                      bias_se = round(bias_se, ndp_bias),
                      sd = round(sd_est, ndp),
                      power = power,
                      power_se = round(power_se, ndp_se),
                      coverage = coverage,
                      coverage_se = round(coverage_se, ndp_se),
                      type1 = type1,
                      mn_ci_width = round(mn_ci_width, ndp))
print(tab_res)



dimnames(sim_curr$E_Y)

# pval_df
# print(colMeans(pval_df < .05))
# par(mfrow = c(1, 2))
# plot(pval_df$fac_0_longit_0_pval, pval_df$fac_0_longit_1_pval)
# abline(0, 1)
# 
# plot(pval_df$fac_0_longit_0_pval, pval_df$fac_1_longit_1_pval)
# abline(0, 1)
#



 





















# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
#       mn_fixed <- summ_fac_1_longit_1$coefficients$fixed
#       cov_fixed <- summ_fac_1_longit_1$varFix
#       se_fixed <- sqrt(diag(cov_fixed))
#       
#       for (arm_curr in arm_all) {
#         if (estimate_type == "contrasts_vs_placebo") {
#           if (arm_curr == placebo_label) {
#             coef_names <- paste0("time_fac", t_unique)
#           } else {
#             coef_names <- paste0("time_fac", t_unique, ":arm", arm_curr)
#           }
#         }
#         if (estimate_type == "treatment_means") {
#           coef_names <- paste0("arm", arm_curr, ":time_fac", t_unique)
#         }
#         m_coef_mat[, arm_curr] <- mn_fixed[coef_names]
#         u_coef_mat[, arm_curr] <-mn_fixed[coef_names] + 2 * se_fixed[coef_names]
#         l_coef_mat[, arm_curr] <- mn_fixed[coef_names] - 2 * se_fixed[coef_names]
#       }
#       res <- list(m = m_coef_mat,
#                       u = u_coef_mat,
#                       l = l_coef_mat)
#     # }
# 
# 
# 
#   # }
#   analysis_results_list[[analysis_curr]] <- results_list
#   analysis_results_list[[analysis_curr]]$mfd_object <- mfd_object_list[[analysis_curr]]
#   analysis_results_list[[analysis_curr]]$arm_all <- arm_all
#   analysis_results_list[[analysis_curr]]$arm_no_plac <- arm_no_plac
#   
# # }
# 

# saveRDS(analysis_results_list, file = file.path("output", "analysis_results_list.RDS"))


# X <- model.matrix(fixed_form_y, model.frame(fixed_form_y, data = results_df))
# ncol(X)
# qr(X)$rank




