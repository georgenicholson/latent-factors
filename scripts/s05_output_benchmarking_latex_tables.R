#######################################
# Source functions
fn_files_to_source <- list.files("functions/", full.names = TRUE)
for (file_curr in fn_files_to_source) {
  source(file_curr)
}

control <- get_control_parameters_simulations(parameter_set = 1)


par_sets_include_in_tab <- names(control$parameter_change_list)

######################################################
# Create list of tables and process using x

latex_for_problem <- "\\textcolor{red}"
latex_for_best <- "\\textbf"
res_list <- list()
new_script <- TRUE
sims_include <- c(1, 2, 5)
for(i in 1:length(control$benchmarking_table_split)) {
  sims_include <- control$benchmarking_table_split[[i]]
  benchmarking_results_file <- control$benchmarking_results_files[i]
  for (j in sims_include) {
      resl <- readRDS(file = file.path("output", 
                paste0("resl_includeVAR_", TRUE, "benchmarking_parameter_set_", j, ".RDS")))
      str(resl$non_null, m = 1)
      fabian_results_dir <- "data/python_deep_learning_output"
      for (meth in c("RNN", "LSTM", "GRU")) {
        resl$non_null[[meth]] <- list()
        resl$non_null[[meth]]$y <- list()
        file_text <- switch(meth, 
                            RNN = "naive",
                            LSTM = "lstm",
                            GRU = "gru")
        resl$non_null[[meth]]$y$est <- 
          read.csv(file.path(fabian_results_dir, paste0("zpkqra53_", file_text, "_", j, ".csv")))$estimated_treatment_effect     
        resl$non_null[[meth]]$run_time <- 
          read.csv(file.path(fabian_results_dir, paste0("zpkqra53_", file_text, "_run_times.csv")))$run_time     
      }

      mae <- sapply(resl[["non_null"]][control$model_names], function(x) mean(abs(x$y$est - resl[["non_null"]]$true_estimand), 
                    na.rm = T))
      mse <- sapply(resl[["non_null"]][control$model_names], function(x) mean((x$y$est - resl[["non_null"]]$true_estimand)^2, 
                    na.rm = T))
      mse_se <- sapply(resl[["non_null"]][control$model_names], function(x) sd((x$y$est - resl[["non_null"]]$true_estimand)^2, 
                                                                               na.rm = T)  / sqrt(sum(!is.na(x$y$est))))
      bias_est <- sapply(resl[["non_null"]][control$model_names], function(x) mean(x$y$est - resl[["non_null"]]$true_estimand, 
                                                                                   na.rm = T))
      bias_v <- sapply(resl[["non_null"]][control$model_names], function(x) x$y$est - resl[["non_null"]]$true_estimand)
      
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
                            mse = format_fun(mse, ndp_precise),
                            mse_se = format_fun(mse_se, ndp),
                            rmse = format_fun(sqrt(mse), ndp_precise),
                            bias = format_fun(bias_est, ndp),
                            bias_se = format_fun(bias_se, ndp),
                            sd = format_fun(sd_est, ndp_precise),
                            sd_se = format_fun(sd_se, ndp_precise),
                            power = format_fun(power, ndp_precise),
                            power_se = format_fun(power_se, ndp),
                            type1 = format_fun(type1, ndp),
                            type1_se = format_fun(type1_se, ndp),
                            bias_problem = abs(bias_est) > 1.96 * bias_se,
                            type1_problem = abs(type1 - control$pval_thresh) > 1.96 * type1_se,
                            coverage_problem = abs(coverage - (1 - control$pval_thresh)) > 1.96 * coverage_se)
      tab_res$rmse[which.min(sqrt(mse))] <- paste0(latex_for_best, "{", tab_res$rmse[which.min(sqrt(mse))], "}")
      tab_res$mse_in <- tab_res$mse
      tab_res$mse[which.min(mse)] <- paste0(latex_for_best, "{", tab_res$mse[which.min(mse)], "}")
      res_list[[j]] <- tab_res
  }
  
  
  
  
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
  nam_out <- c('Method', 'Input', 'MV', 'Longit', bias_nam, sd_nam, mse_nam)
  
  res_list_out <- lapply(res_list[sims_include], function(x) {
                                x$Method <- control$methods;
                                x$Input <- paste0("\\dataset{", control$input, "}");
                                x$Longit <- ifelse(control$input %in% 3:4, "\\checkmark", "");
                                x$MV <- ifelse(control$input %in% c(2, 4), "\\checkmark", "");
                                x[mse_nam] <- paste0(x$mse, se_fun(x$mse_se));
                                x[bias_nam] <- paste0(ifelse(x$bias_problem, latex_for_problem, ""),
                                                          "{", x$bias, se_fun(x$bias_se), "}");
                                x[sd_nam] <- paste0(x$sd, se_fun(x$sd_se));
                                x[pow_nam] <- paste0(x$power, se_fun(x$power_se));
                                x[coverage_nam] <- paste0(ifelse(x$coverage_problem, latex_for_problem, ""),
                                                            "{", x$coverage, se_fun(x$coverage_se), "}");
                                x[type1_nam] <- paste0(ifelse(x$type1_problem, latex_for_problem, ""),
                                                          "{", x$type1, se_fun(x$type1_se), "}");
                                return(x[, nam_out])})
  
  subtable_headings <- paste0(names(control$parameter_change_list)[sims_include], " = ", 
         sapply(control$parameter_change_list[sims_include], function(x) x[1]))
  subtable_headings[subtable_headings == "$\\Theta_1$: Default parameters = NULL"] <- "$\\Theta_1$: Default parameters (see Table~\\ref{tab:parameter_defaults})"
  subtable_headings <- paste0("(", letters[1:length(control$parameter_change_list[sims_include])], ") ", subtable_headings)
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
  cat(sim_results_latex, file = benchmarking_results_file)
}

#######################
# Output runtime table

j <- 1
resl <- readRDS(file = file.path("output", 
                                 paste0("resl_includeVAR_", TRUE, "benchmarking_parameter_set_", j, ".RDS")))

  
fabian_results_dir <- "sim_results_from_Fabian"
for (meth in c("RNN", "LSTM", "GRU")) {
  resl$non_null[[meth]] <- list()
  resl$non_null[[meth]]$y <- list()
  file_text <- switch(meth, 
                      RNN = "naive",
                      LSTM = "lstm",
                      GRU = "gru")
    resl$non_null[[meth]]$run_time <- 
      read.csv(file.path(fabian_results_dir, paste0("zpkqra53_", file_text, "_run_times.csv")))$run_time     
}


t_df <- data.frame(Method = control$methods,
                   Input = paste0("\\dataset{", control$input, "}"),
                   Longit = ifelse(control$input %in% 3:4, "\\checkmark", ""),
                   MV = ifelse(control$input %in% c(2, 4), "\\checkmark", ""),
                   Hardware = c(rep("Intel Xeon", length(control$input) - 3),
                                rep("NVIDIA A100", 3)),
                   'Complexity' = control$complexity,
                  'Run time (s)' = sapply(resl$non_null[control$model_names], function(x) 
                          round(mean(x$run_time, na.rm = T), 1))[control$model_names],
                   check.names = F) 
tabout <- xtable::xtable(t_df, digits = 1)
complexity_out <- print(tabout, 
                        table.placement = "hp", 
                        colnames.format = "multiple", 
                        caption.placement = "top", 
                        sanitize.text.function = function(x){x}, 
                        include.rownames = F, 
                        floating = FALSE)


dir.create("tables", showWarnings = FALSE)
cat(complexity_out, file = control$benchmarking_runtime_file)



######################################################
# Export numbers for text

# MSE comp sim table (a) 
runt <- t_df$`Run time (s)`
names(runt) <- t_df$Method
lm_var_fold <- formatC(runt["VAR"] / runt["\\lm{4}"], digits = 0, format = "f")
lm_rnn_fold <- formatC(mean(c(runt["RNN"],  runt["LSTM"],  runt["GRU"])) / runt["\\lm{4}"], digits = 0, format = "f")
lm4_secs <- formatC(runt["\\lm{4}"], digits = 1, format = "f")
gru_secs <- formatC(runt["GRU"], digits = 1, format = "f")
lstm_secs <- formatC(runt["LSTM"], digits = 1, format = "f")
rnn_secs <- formatC(runt["RNN"], digits = 1, format = "f")
var_secs <- formatC(runt["VAR"], digits = 1, format = "f")
save.num <- c("lm_var_fold", "lm_rnn_fold", "lm4_secs", "gru_secs", "lstm_secs", "rnn_secs", "var_secs")
for (numc in save.num)
  write.table(eval(as.name(numc)), 
              file = file.path(control$numbers_out_dir, paste0(numc, ".txt")),
              col.names = F, row.names = F, quote = F)



