#######################################
# Source functions
fn_files_to_source <- list.files("functions/", full.names = TRUE)
for (file_curr in fn_files_to_source) {
  source(file_curr)
}

control <- get_control_parameters_simulations(parameter_set = 1)

par_sets_include_in_tab <- c("Default parameters", "Increased magnitude of (negative) effect size $\\mu_{\\text{active,16}}$", 
  "Decreased effect size", "Increased $\\rhoARNok$", "Decreased $\\rhoARNok$", 
  "Increased $\\sigmaARNok$", "Increased $\\sigmaRINok$", "Increased $\\sigmaY$")[-3]


######################################################
# Create list of tables and process using x

latex_for_problem <- "\\textcolor{red}"
latex_for_best <- "\\textbf"

# model_names <- c("fac_0_longit_0", "fac_1_longit_0", "fac_0_longit_1", "fac_1_longit_1")
res_list <- list()
new_script <- TRUE
j <- 1#for (j in 1:control$n_parameter_sets) {
sims_include <- c(1, 2, 5)
for (j in sims_include) {
  # if (new_script) {
    resl <- readRDS(file = file.path("output", 
                                     paste0("resl_benchmarking_parameter_set_", j, ".RDS")))
    str(resl$non_null, m = 1)
    mae <- sapply(resl[["non_null"]][control$model_names], function(x) mean(abs(x$y$est - resl[["non_null"]]$true_estimand), 
                  na.rm = T))
    mse <- sapply(resl[["non_null"]][control$model_names], function(x) mean((x$y$est - resl[["non_null"]]$true_estimand)^2, 
                  na.rm = T))
    mse_se <- sapply(resl[["non_null"]][control$model_names], function(x) sd((x$y$est - resl[["non_null"]]$true_estimand)^2, 
                                                                             na.rm = T)  / sqrt(sum(!is.na(x$y$est))))
    bias_est <- sapply(resl[["non_null"]][control$model_names], function(x) mean(x$y$est - resl[["non_null"]]$true_estimand, 
                                                                                 na.rm = T))
    bias_se <- sapply(resl[["non_null"]][control$model_names], function(x) sd(x$y$est - resl[["non_null"]]$true_estimand, 
                                                                              na.rm = T) / sqrt(sum(!is.na(x$y$est))))
    sd_est <- sapply(resl[["non_null"]][control$model_names], function(x) sd(x$y$est - resl[["non_null"]]$true_estimand, 
                                                                             na.rm = T))
    sd_se <- sapply(resl[["non_null"]][control$model_names], function(x) sd((x$y$est - resl[["non_null"]]$true_estimand)^2, 
                                                                             na.rm = T) / sqrt(control$n_its - 1))
    power <- sapply(resl[["non_null"]][control$model_names], function(x) mean(x$y$pval < control$pval_thresh, 
                                                                              na.rm = T))
    power_se <- sapply(resl[["non_null"]][control$model_names], function(x) sd(x$y$pval < control$pval_thresh, 
                                                                               na.rm = T) / sqrt(sum(!is.na(x$y$est))))
    type1 <- sapply(resl[["null"]][control$model_names], function(x) mean(x$y$pval < control$pval_thresh, 
                                                                          na.rm = T))
    type1_se <- sapply(resl[["null"]][control$model_names], function(x) sd(x$y$pval < control$pval_thresh, 
                                                                           na.rm = T) / sqrt(sum(!is.na(x$y$est))))
    coverage <- sapply(resl[["non_null"]][control$model_names], function(x) mean(x$y$low <= resl[["non_null"]]$true_estimand & 
                                                                           x$y$upp >= resl[["non_null"]]$true_estimand, 
                                                                           na.rm = T))
    coverage_se <- sapply(resl[["non_null"]][control$model_names], function(x) sd(x$y$low <= resl[["non_null"]]$true_estimand & 
                                                                           x$y$upp >= resl[["non_null"]]$true_estimand, 
                                                                           na.rm = T) / sqrt(sum(!is.na(x$y$est))))
    mn_ci_width <- sapply(resl[["non_null"]][control$model_names], function(x) mean(x$y$upp - x$y$low, 
                                                                                    na.rm = T))
    mn_ci_width_se <- sapply(resl[["non_null"]][control$model_names], function(x) sd(x$y$upp - x$y$low, 
                                                                                     na.rm = T) / sqrt(sum(!is.na(x$y$est))))
    
    ndp <- 2
    ndp_precise <- 2
    format_fun <- function(x, nsf = 3) {
      sprintf(paste0('%#.', nsf, 'g'), x)
    }
    tab_res <- data.frame(meth = names(mse),
                          input = control$input,
                          # fac = c(F, T, F, T),
                          # longit = c(F, F, T, T),
                          # model = c("M1", "M2", "M3", "M4"),
                          # mae = format_fun(mae, ndp_precise),
                          mse = format_fun(mse, ndp_precise),
                          mse_se = format_fun(mse_se, ndp),
                          rmse = format_fun(sqrt(mse), ndp_precise),
                          bias = format_fun(bias_est, ndp),
                          # bias2 = format_fun(bias_est^2, ndp),
                          bias_se = format_fun(bias_se, ndp),
                          sd = format_fun(sd_est, ndp_precise),
                          sd_se = format_fun(sd_se, ndp_precise),
                          power = format_fun(power, ndp_precise),
                          power_se = format_fun(power_se, ndp),
                          # coverage = format_fun(coverage, ndp_precise),
                          # coverage_se = format_fun(coverage_se, ndp),
                          type1 = format_fun(type1, ndp),
                          type1_se = format_fun(type1_se, ndp),
                          # mn_ci_width = format_fun(mn_ci_width, ndp_precise),
                          # mn_ci_width_se = format_fun(mn_ci_width_se, ndp),
                          bias_problem = abs(bias_est) > 1.96 * bias_se,
                          type1_problem = abs(type1 - control$pval_thresh) > 1.96 * type1_se,
                          coverage_problem = abs(coverage - (1 - control$pval_thresh)) > 1.96 * coverage_se)
    tab_res$rmse[which.min(sqrt(mse))] <- paste0(latex_for_best, "{", tab_res$rmse[which.min(sqrt(mse))], "}")
    tab_res$mse_in <- tab_res$mse
    tab_res$mse[which.min(mse)] <- paste0(latex_for_best, "{", tab_res$mse[which.min(mse)], "}")
    tab_res$bias[which.min(bias_est^2)] <- paste0(latex_for_best, "{", tab_res$bias[which.min(bias_est^2)], "}")
    tab_res$sd[which.min(tab_res$sd)] <- paste0(latex_for_best, "{", tab_res$sd[which.min(tab_res$sd)], "}")
    # tab_res$power[which.max(power)] <- paste0(latex_for_best, "{", tab_res$power[which.max(power)], "}")
    # tab_res$mn_ci_width[which.min(mn_ci_width)] <- paste0(latex_for_best, "{", tab_res$mn_ci_width[which.min(mn_ci_width)], "}")
    res_list[[j]] <- tab_res
}



  # } else {
  #   res_list[[j]] <- readRDS(file = file.path("output", paste0("tab_res_parameter_set_", j, ".RDS")))
  # }
# }
# names(res_list) <- names(control$parameter_change_list)

se_fun <- function(x) {
  paste0(" $_{(", x, ")}$")
}
se_fun_title <- function(x) {
  paste0(" $_{(", x, ")}$")
}

pow_nam <- paste0('Power', se_fun_title("\\text{SE}"))
mse_nam <- paste0('MSE', se_fun_title("\\text{SE}"))
sd_nam <- paste0('SD', se_fun_title("\\text{SE}"))
mn_ci_width_nam <- paste0('CI width', se_fun_title("\\text{SE}"))
coverage_nam <- paste0('Coverage', se_fun_title("\\text{SE}"))
bias_nam <- paste0('Bias', se_fun_title("\\text{SE}"))
type1_nam <- paste0('Type I error', se_fun_title("\\text{SE}"))
# nam_out <- c('Model', bias_nam, coverage_nam, mse_nam, mn_ci_width_nam, pow_nam)
nam_out <- c('Method', 'Input', 'MV', 'Longit', mse_nam, bias_nam, sd_nam)#, type1_nam, pow_nam)
x <- par_sets_include_in_tab[2]

res_list_out <- lapply(res_list[sims_include], function(x) {
                              x$Method <- control$methods;
                              x$Input <- paste0("\\dataset{", control$input, "}");
                              x$Longit <- ifelse(control$input %in% 3:4, "\\checkmark", "");
                              x$MV <- ifelse(control$input %in% c(2, 4), "\\checkmark", "");
                              x[mse_nam] <- paste0(x$mse, se_fun(x$mse_se));
                              x[bias_nam] <- paste0(ifelse(x$bias_problem, latex_for_problem, ""), 
                                                      "{", x$bias, se_fun(x$bias_se), "}");
                              x[sd_nam] <- paste0(x$sd, se_fun(x$sd_se));
                              # x$`Power $\\uparrow$` <- x$power;
                              x[pow_nam] <- paste0(x$power, se_fun(x$power_se));
                              # x$`Coverage` <- x$coverage;
                              x[coverage_nam] <- paste0(ifelse(x$coverage_problem, latex_for_problem, ""),
                                                          "{", x$coverage, se_fun(x$coverage_se), "}");
                              x[type1_nam] <- paste0(ifelse(x$type1_problem, latex_for_problem, ""),
                                                        "{", x$type1, se_fun(x$type1_se), "}");
                              # x$SD <- paste0(x$sd);
                              # x$`MSE $\\downarrow$` <- x$mse;
                              # x$`Type I err` <- x$type1;
                              # x$`Type I err (SE)` <- paste0(ifelse(x$type1_problem, latex_for_problem, ""), 
                              #                               "{", x$type1, "$_{(", x$type1_se, ")}$", "}");
                              # x[mn_ci_width_nam] <- paste0(x$mn_ci_width, se_fun(x$mn_ci_width_se));
                              return(x[, nam_out])})

subtable_headings <- paste0(names(control$parameter_change_list[par_sets_include_in_tab]), " = ", 
       sapply(control$parameter_change_list[par_sets_include_in_tab], function(x) x[1]))
subtable_headings[subtable_headings == "Default parameters = NULL"] <- "Default parameters (see Table~\\ref{tab:parameter_defaults})"
subtable_headings <- paste0("(", letters[1:length(par_sets_include_in_tab)], ") ", subtable_headings)
attr(res_list_out, "subheadings") <- subtable_headings


if (length(res_list_out) == 1) {
  tabout <- xtable::xtable(res_list_out[[1]])
  sim_results_latex <- print(tabout, 
                             table.placement = "hp", 
                             colnames.format = "multiple", 
                             caption.placement = "top", 
                             sanitize.text.function = function(x){x}, 
                             include.rownames = F, 
                             floating = FALSE)
} else {
  xList <- xtable::xtableList(res_list_out)
  sim_results_latex <- print(xList,
                             table.placement = "hp",
                             colnames.format = "multiple",
                             caption.placement = "top",
                             sanitize.text.function = function(x){x},
                             include.rownames = F,
                             floating = FALSE)
  
}


dir.create("tables", showWarnings = FALSE)
cat(sim_results_latex, file = control$benchmarking_results_file)



# Paste below in terminal to copy to Overleaf
# tar xvzf '/mnt/c/Users/nicho/Dropbox/GitHub_Dropbox/latent-factors/zipped_tables.tar.gz' -C '/mnt/c/Users/nicho/Dropbox/Apps/Overleaf/LatentFactors_OxfordNovartis (J Biomedical Informatics)/Figures/BDI_plot_export' 

# 
# browseURL(tab_export_files_zipped) # this downloads locally to Downloads folder on a PC
# 
# 
# 
# # 
# # 
# 
# 
# 
# dir.create("text_numbers", showWarnings = FALSE)
# numbers_out <- "text_numbers/"
# 
# ######################################################
# # Export numbers for text
# 
# # MSE comp sim table (a) 
# tab_curr <- res_list$`Default parameters`
# m2m4_mse_a <- paste0("MSE$_\\text{M2}=", tab_curr$mse_in[tab_curr$model == "M2"], "$ and ",  
#                      "MSE$_\\text{M4}=", tab_curr$mse_in[tab_curr$model == "M4"], "$")
# m1m3_mse_a <- paste0("MSE$_\\text{M1}=", tab_curr$mse_in[tab_curr$model == "M1"], "$ and ",  
#                      "MSE$_\\text{M3}=", tab_curr$mse_in[tab_curr$model == "M3"], "$")
# save.num <- c("m2m4_mse_a", "m1m3_mse_a")
# for (numc in save.num)
#   write.table(eval(as.name(numc)), 
#               file = file.path("text_numbers", paste0(numc, ".txt")),
#               col.names = F, row.names = F, quote = F)
# 
# # MSE % change from a to e
# mse_a <- res_list$`Default parameters`$mse_in
# mse_e <- res_list$`Increased $\\sigmaARNok$`$mse_in
# mse_a_to_e <- paste(round(range(as.numeric(mse_e) / as.numeric(mse_a)), 2), collapse = " to ")
# for (numc in "mse_a_to_e")
#   write.table(eval(as.name(numc)), 
#               file = file.path("text_numbers", paste0(numc, ".txt")),
#               col.names = F, row.names = F, quote = F)
# 
# # MSE change from a to f for longit vs non longit models
# tab_a <- res_list$`Default parameters`
# tab_f <- res_list$`Increased $\\sigmaRINok$`
# mse_a_longit <- as.numeric(tab_a$mse_in[tab_a$longit])
# mse_f_longit <- as.numeric(tab_f$mse_in[tab_f$longit])
# mse_a_non_longit <- as.numeric(tab_a$mse_in[!tab_a$longit])
# mse_f_non_longit <- as.numeric(tab_f$mse_in[!tab_f$longit])
# mse_a_to_f_longit <- paste(round(range(as.numeric(mse_f_longit) / as.numeric(mse_a_longit)), 2), collapse = " and ")
# mse_a_to_f_non_longit <- paste(round(range(as.numeric(mse_f_non_longit) / as.numeric(mse_a_non_longit)), 2), collapse = " and ")
# for (numc in c("mse_a_to_f_longit", "mse_a_to_f_non_longit"))
#   write.table(eval(as.name(numc)), 
#               file = file.path("text_numbers", paste0(numc, ".txt")),
#               col.names = F, row.names = F, quote = F)
# 
# 
# # MSE change from a to g for longit vs non longit models
# tab_a <- res_list$`Default parameters`
# tab_g <- res_list$`Increased $\\sigmaY$`
# mse_a_fac <- as.numeric(tab_a$mse_in[tab_a$fac])
# mse_e_fac <- as.numeric(tab_g$mse_in[tab_g$fac])
# mse_a_non_fac <- as.numeric(tab_a$mse_in[!tab_a$fac])
# mse_e_non_fac <- as.numeric(tab_g$mse_in[!tab_g$fac])
# mse_e <- res_list$`Increased $\\sigmaARNok$`$mse_in
# mse_a_to_g_fac <- paste(round(range(as.numeric(mse_e_fac) / as.numeric(mse_a_fac)), 2), collapse = " and ")
# mse_a_to_g_non_fac <- paste(round(range(as.numeric(mse_e_non_fac) / as.numeric(mse_a_non_fac)), 2), collapse = " and ")
# for (numc in c("mse_a_to_g_fac", "mse_a_to_g_non_fac"))
#   write.table(eval(as.name(numc)), 
#               file = file.path("text_numbers", paste0(numc, ".txt")),
#               col.names = F, row.names = F, quote = F)
# 
# 
# 
# ######################################################
# # Export the full set of tables and emax plot in tar.gz format
# source("scripts/emax_curves_demo_plot.R")
# tab_export_files <- c(default_sim_table_file, sim_results_file, emax_example_plot_file, numbers_out)
# tab_export_files_zipped <- zip::zip(zipfile = file.path("zipped_tables.zip"),
#                                     files = tab_export_files,
#                                     mode = "cherry-pick")
# tab_export_files_zipped <- file.path("zipped_tables.tar.gz")
# tar(tarfile = tab_export_files_zipped,
#     files = tab_export_files,
#     compression = 'gzip',
#     tar = "tar")
# browseURL(tab_export_files_zipped) # this downloads locally to Downloads folder on a PC
# 
# 

