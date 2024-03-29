#######################################
# Source functions
fn_files_to_source <- list.files("functions/", full.names = TRUE)
for (file_curr in fn_files_to_source) {
  source(file_curr)
}

control <- get_control_parameters_simulations(parameter_set = 1)

par_sets_include_in_tab <- c("$\\Theta_1$: Default parameters", 
                             "$\\Theta_2$: Increased magnitude of (negative) effect size $\\mu_{\\text{active,16}}$", 
                             "$\\Theta_3$: Increased $\\rhoARNok$", 
                             "$\\Theta_4$: Decreased $\\rhoARNok$", 
                             "$\\Theta_5$: Increased $\\sigmaARNok$", 
                             "$\\Theta_6$: Increased $\\sigmaRINok$", 
                             "$\\Theta_7$: Increased $\\sigmaY$")


######################################################
# Create parameter table
tab_out <- data.frame()
tab_out <- rbind(tab_out, data.frame(var_name = "n_subj", Parameter = "$N$", Description = "Total number of trial participants", round = 0,
                                     Defined = "Section~\\ref{sec:notation}"))
tab_out <- rbind(tab_out, data.frame(var_name = "n_var", Parameter = "$P$", Description = "Dimension of $\\B{Y}$, number of measurements", round = 0,
                                     Defined = "Section~\\ref{sec:notation}"))
tab_out <- rbind(tab_out, data.frame(var_name = "n_fac", Parameter = "$K$", Description = "Dimension of $\\B{Z}$, number of latent variables", round = 0,
                                     Defined = "Section~\\ref{sec:notation}"))
tab_out <- rbind(tab_out, data.frame(var_name = "time_of_half_max_effect", Parameter = "$EC_{50}$", Description = "Time of half-max response", round = 0,
                                     Defined = "Eqn~\\eqref{eq:emax_formula}"))
tab_out <- rbind(tab_out, data.frame(var_name = "size_of_final_week_effect_placebo", 
                                     Parameter = "$\\mu_{\\text{placebo},16}$", 
                                     Description = "Placebo arm mean response, week 16", round = 1,
                                     Defined = "Eqn~\\eqref{eq:emax_formula}"))
tab_out <- rbind(tab_out, data.frame(var_name = "size_of_final_week_effect_active", 
                                     Parameter = "$\\mu_{\\text{active},16}$", 
                                     Description = "Active arm mean response, week 16", round = 1,
                                     Defined = "Eqn~\\eqref{eq:emax_formula}"))
tab_out <- rbind(tab_out, data.frame(var_name = "theta_RI_sd", Parameter = "$\\sigmaRINok$", Description = "Random intercept SD for $\\BZ$", 
                                     round = 2,
                                     Defined = "Eqn~\\eqref{eq:Sigdef1}"))
tab_out <- rbind(tab_out, data.frame(var_name = "theta_RI_sd", Parameter = "$\\sigmaRINok$", Description = "Random intercept SD for $\\BZ$", 
                                     round = 2,
                                     Defined = "Eqn~\\eqref{eq:Sigdef1}"))
tab_out <- rbind(tab_out, data.frame(var_name = "theta_resid_sd", Parameter = "$\\sigmaARNok$", Description = "AR1 SD for $\\BZ$", round = 2,
                                     Defined = "Eqn~\\eqref{eq:Sigdef1}"))
tab_out <- rbind(tab_out, data.frame(var_name = "cor_ar1_lag", Parameter = "$\\rhoARNok$", Description = "AR1 correlation for $\\BZ$", round = 2,
                                     Defined = "Eqn~\\eqref{eq:Sigdef1}"))
tab_out <- rbind(tab_out, data.frame(var_name = "y_resid_sd", Parameter = "$\\sigmaY$", Description = "Residual SD for $\\B{Y}$", round = 2,
                                     Defined = "Eqn~\\eqref{eq:simple_stage1}"))
tab_out <- rbind(tab_out, data.frame(var_name = "A", Parameter = "$\\B{A}$", Description = "Sparse $P\\times K$ loadings matrix", round = 2,
                                     Defined = "Eqn~\\eqref{eq:simple_stage1}"))

tab_out$`Default value` <- NA
for (j in which(tab_out$var_name != "A")) {
  tab_out$`Default value`[j] <- formatC(control[[tab_out$var_name[j]]], digits = tab_out$round[j], format = "f")
}
tab_out
tab_out[match("A", tab_out$var_name), "Default value"] <- "Figure~\\ref{fig:loadings_comparison}(a)"


######################################################
# Take subset of columns and create txt file with latex code
tab_print <- tab_out[, c("Parameter", "Description", "Defined", "Default value")]
xtab_tab_print <- print(xtable::xtable(tab_print, 
             label = "tab:parameter_defaults", 
             align = rep("l", ncol(tab_print) + 1),
             caption = ""),
      caption.placement = "top", 
      sanitize.text.function = function(x){x},
      include.rownames = F, 
      floating = FALSE)
dir.create("tables", showWarnings = FALSE)
cat(xtab_tab_print, file = control$default_sim_table_file)


######################################################
# Create list of tables and process using x

latex_for_problem <- "\\textcolor{red}"
latex_for_best <- "\\textbf"

model_names <- c("fac_0_longit_0", "fac_1_longit_0", "fac_0_longit_1", "fac_1_longit_1")
res_list <- list()
new_script <- TRUE
for (j in 1:control$n_parameter_sets) {
  if (new_script) {
    resl <- readRDS(file = file.path("output", 
                                     paste0("resl_parameter_set_", j, ".RDS")))
    mae <- sapply(resl[["non_null"]][model_names], function(x) mean(abs(x$y$est - resl[["non_null"]]$true_estimand)))
    mse <- sapply(resl[["non_null"]][model_names], function(x) mean((x$y$est - resl[["non_null"]]$true_estimand)^2))
    rmse <- sqrt(mse)
    mse_se <- sapply(resl[["non_null"]][model_names], function(x) sd((x$y$est - resl[["non_null"]]$true_estimand)^2)  / sqrt(control$n_its))
    bias_est <- sapply(resl[["non_null"]][model_names], function(x) mean(x$y$est - resl[["non_null"]]$true_estimand))
    bias_v <- sapply(resl[["non_null"]][model_names], function(x) x$y$est - resl[["non_null"]]$true_estimand)
    bias_se <- sapply(resl[["non_null"]][model_names], function(x) sd(x$y$est - resl[["non_null"]]$true_estimand) / sqrt(control$n_its))
    sd_est <- sapply(resl[["non_null"]][model_names], function(x) sd(x$y$est - resl[["non_null"]]$true_estimand))
    power <- sapply(resl[["non_null"]][model_names], function(x) mean(x$y$pval < control$pval_thresh))
    power_se <- sapply(resl[["non_null"]][model_names], function(x) sd(x$y$pval < control$pval_thresh) / sqrt(control$n_its))
    type1 <- sapply(resl[["null"]][model_names], function(x) mean(x$y$pval < control$pval_thresh))
    type1_se <- sapply(resl[["null"]][model_names], function(x) sd(x$y$pval < control$pval_thresh) / sqrt(control$n_its))
    coverage <- sapply(resl[["non_null"]][model_names], function(x) mean(x$y$low <= resl[["non_null"]]$true_estimand & 
                                                                           x$y$upp >= resl[["non_null"]]$true_estimand))
    coverage_se <- sapply(resl[["non_null"]][model_names], function(x) sd(x$y$low <= resl[["non_null"]]$true_estimand & 
                                                                           x$y$upp >= resl[["non_null"]]$true_estimand) / sqrt(control$n_its))
    mn_ci_width <- sapply(resl[["non_null"]][model_names], function(x) mean(x$y$upp - x$y$low))
    mn_ci_width_se <- sapply(resl[["non_null"]][model_names], function(x) sd(x$y$upp - x$y$low) / sqrt(control$n_its))
    
    ndp <- 2
    ndp_precise <- 2
    format_fun <- function(x, nsf = 3) {
      sprintf(paste0('%#.', nsf, 'g'), x)
    }
    tab_res <- data.frame(meth = names(mse),
                          fac = c(F, T, F, T),
                          longit = c(F, F, T, T),
                          model = c("M1", "M2", "M3", "M4"),
                          input = 1:4,
                          mae = format_fun(mae, ndp_precise),
                          mse = format_fun(mse, ndp_precise),
                          mse_se = format_fun(mse_se, ndp),
                          rmse = format_fun(sqrt(mse), ndp_precise),
                          bias = format_fun(bias_est, ndp),
                          bias_se = format_fun(bias_se, ndp),
                          sd = format_fun(sd_est, ndp_precise),
                          power = format_fun(power, ndp_precise),
                          power_se = format_fun(power_se, ndp),
                          coverage = format_fun(coverage, ndp_precise),
                          coverage_se = format_fun(coverage_se, ndp),
                          type1 = format_fun(type1, ndp),
                          type1_se = format_fun(type1_se, ndp),
                          mn_ci_width = format_fun(mn_ci_width, ndp_precise),
                          mn_ci_width_se = format_fun(mn_ci_width_se, ndp),
                          bias_problem = abs(bias_est) > 1.96 * bias_se,
                          type1_problem = type1 - control$pval_thresh > 1.96 * type1_se,
                          coverage_problem = (1 - control$pval_thresh) - coverage > 1.96 * coverage_se)
    tab_res$rmse[which.min(sqrt(mse))] <- paste0(latex_for_best, "{", rmse[which.min(sqrt(mse))], "}")
    tab_res$mse_in <- tab_res$mse
    tab_res$mse[which.min(mse)] <- paste0(latex_for_best, "{", tab_res$mse[which.min(mse)], "}")
    tab_res$power[which.max(power)] <- paste0(latex_for_best, "{", tab_res$power[which.max(power)], "}")
    tab_res$mn_ci_width[which.min(mn_ci_width)] <- paste0(latex_for_best, "{", tab_res$mn_ci_width[which.min(mn_ci_width)], "}")
    res_list[[j]] <- tab_res
    
  } else {
    res_list[[j]] <- readRDS(file = file.path("output", paste0("tab_res_parameter_set_", j, ".RDS")))
  }
}
names(res_list) <- names(control$parameter_change_list)

se_fun <- function(x) {
  paste0(" $_{(", x, ")}$")
}
se_fun_title <- function(x) {
  paste0(" $_{(", x, ")}$")
}

pow_nam <- paste0('Power', se_fun_title("\\text{SE}"))
mse_nam <- paste0('MSE', se_fun_title("\\text{SE}"))
mn_ci_width_nam <- paste0('CI width', se_fun_title("\\text{SE}"))
coverage_nam <- paste0('Coverage', se_fun_title("\\text{SE}"))
bias_nam <- paste0('Bias', se_fun_title("\\text{SE}"))
type1_nam <- paste0('Type I error', se_fun_title("\\text{SE}"))

nam_out <- c('Method', 'Input', 'MV', 'Longit', bias_nam, type1_nam, coverage_nam, pow_nam, mse_nam)


res_list_out <- lapply(res_list[par_sets_include_in_tab], function(x) {
                              x$Model <- paste0("M", 1:4);
                              x$Method <- paste0("\\lm{", 1:4, "}");
                              x$Input <- paste0("\\dataset{", 1:4, "}");
                              x$Longit <- c("", "", "\\checkmark", "\\checkmark");
                              x$MV <- c("", "\\checkmark", "", "\\checkmark");
                              x[bias_nam] <- paste0(ifelse(x$bias_problem, latex_for_problem, ""), 
                                                      "{", x$bias, se_fun(x$bias_se), "}");
                              x[mse_nam] <- paste0(x$mse, se_fun(x$mse_se));
                              x$`Power $\\uparrow$` <- x$power;
                              x[pow_nam] <- paste0(x$power, se_fun(x$power_se));
                              x$`Coverage` <- x$coverage;
                              x[coverage_nam] <- paste0(ifelse(x$coverage_problem, latex_for_problem, ""), 
                                                          "{", x$coverage, se_fun(x$coverage_se), "}");
                              x[type1_nam] <- paste0(ifelse(x$type1_problem, latex_for_problem, ""),
                                                     "{", x$type1, se_fun(x$type1_se), "}");
                              x$SD <- paste0(x$sd);
                              x$`MSE $\\downarrow$` <- x$mse;
                              x$`Type I err` <- x$type1;
                              x$`Type I err (SE)` <- paste0(ifelse(x$type1_problem, latex_for_problem, ""), 
                                                            "{", x$type1, "$_{(", x$type1_se, ")}$", "}");
                              x[mn_ci_width_nam] <- paste0(x$mn_ci_width, se_fun(x$mn_ci_width_se));
                              return(x[, nam_out])})

subtable_headings <- paste0(names(control$parameter_change_list[par_sets_include_in_tab]), " = ", 
       sapply(control$parameter_change_list[par_sets_include_in_tab], function(x) x[1]))
subtable_headings[subtable_headings == "$\\Theta_1$: Default parameters = NULL"] <- 
  "$\\Theta_1$: Default parameters (see Table~\\ref{tab:parameter_defaults})"
subtable_headings <- paste0("(", letters[1:length(par_sets_include_in_tab)], ") ", subtable_headings)
attr(res_list_out, "subheadings") <- subtable_headings
xList <- xtable::xtableList(res_list_out)
sim_results_latex <- print(xList, 
                     table.placement = "hp", 
                     colnames.format = "multiple", 
                     caption.placement = "top", 
                     sanitize.text.function = function(x){x}, 
                     include.rownames = F, 
                     floating = FALSE)
cat(sim_results_latex, file = control$sim_results_file)



dir.create(control$numbers_out, showWarnings = FALSE)

######################################################
# Export numbers for text

# MSE comp sim table (a) 
tab_curr <- res_list$`$\\Theta_1$: Default parameters`
m2m4_mse_a <- paste0("MSE$_{\\text{\\lm{2}}}=", tab_curr$mse_in[tab_curr$model == "M2"], "$ and ",  
                     "MSE$_{\\text{\\lm{4}}}=", tab_curr$mse_in[tab_curr$model == "M4"], "$")
m1m3_mse_a <- paste0("MSE$_{\\text{\\lm{1}}}=", tab_curr$mse_in[tab_curr$model == "M1"], "$ and ",  
                     "MSE$_{\\text{\\lm{3}}}=", tab_curr$mse_in[tab_curr$model == "M3"], "$")
save.num <- c("m2m4_mse_a", "m1m3_mse_a")
for (numc in save.num)
  write.table(eval(as.name(numc)), 
              file = file.path(control$numbers_out_dir, paste0(numc, ".txt")),
              col.names = F, row.names = F, quote = F)

# MSE % change from a to e
mse_a <- res_list$`$\\Theta_1$: Default parameters`$mse_in
mse_e <- res_list$`$\\Theta_5$: Increased $\\sigmaARNok$`$mse_in
mse_a_to_e <- paste(round(base::range(as.numeric(mse_e) / as.numeric(mse_a)), 2), collapse = " to ")
for (numc in "mse_a_to_e")
  write.table(eval(as.name(numc)), 
              file = file.path(control$numbers_out_dir, paste0(numc, ".txt")),
              col.names = F, row.names = F, quote = F)

# MSE change from a to f for longit vs non longit models
tab_a <- res_list$`$\\Theta_1$: Default parameters`
tab_f <- res_list$`$\\Theta_6$: Increased $\\sigmaRINok$`
mse_a_longit <- as.numeric(tab_a$mse_in[tab_a$longit])
mse_f_longit <- as.numeric(tab_f$mse_in[tab_f$longit])
mse_a_non_longit <- as.numeric(tab_a$mse_in[!tab_a$longit])
mse_f_non_longit <- as.numeric(tab_f$mse_in[!tab_f$longit])
mse_a_to_f_longit <- paste(round(base::range(as.numeric(mse_f_longit) / as.numeric(mse_a_longit)), 2), collapse = " and ")
mse_a_to_f_non_longit <- paste(round(base::range(as.numeric(mse_f_non_longit) / as.numeric(mse_a_non_longit)), 2), collapse = " and ")
for (numc in c("mse_a_to_f_longit", "mse_a_to_f_non_longit"))
  write.table(eval(as.name(numc)), 
              file = file.path(control$numbers_out_dir, paste0(numc, ".txt")),
              col.names = F, row.names = F, quote = F)


# MSE change from a to g for longit vs non longit models
tab_a <- res_list$`$\\Theta_1$: Default parameters`
tab_g <- res_list$`$\\Theta_7$: Increased $\\sigmaY$`
mse_a_fac <- as.numeric(tab_a$mse_in[tab_a$fac])
mse_e_fac <- as.numeric(tab_g$mse_in[tab_g$fac])
mse_a_non_fac <- as.numeric(tab_a$mse_in[!tab_a$fac])
mse_e_non_fac <- as.numeric(tab_g$mse_in[!tab_g$fac])
mse_e <- res_list$`Increased $\\sigmaARNok$`$mse_in
mse_a_to_g_fac <- paste(round(base::range(as.numeric(mse_e_fac) / as.numeric(mse_a_fac)), 2), collapse = " and ")
mse_a_to_g_non_fac <- paste(round(base::range(as.numeric(mse_e_non_fac) / as.numeric(mse_a_non_fac)), 2), collapse = " and ")
for (numc in c("mse_a_to_g_fac", "mse_a_to_g_non_fac"))
  write.table(eval(as.name(numc)), 
              file = file.path(control$numbers_out_dir, paste0(numc, ".txt")),
              col.names = F, row.names = F, quote = F)





