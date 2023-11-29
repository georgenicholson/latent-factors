source("scripts/04_plotting_intro.R")


study_plot <- "F2342"
analysis_to_plot <- "PsA_F2342"
subj_unique <- unique(d_time$USUBJID[d_time$study == study_plot])
n_miss<-n_week <- c()

for (subj_plot in subj_unique) {
  d_time_sub <- d_time[which(d_time$USUBJID == subj_plot & d_time$AVISITN %in% c(0, 2, 4, 8, 12, 16)), ]
  ymat <- d_time_sub[, nms_params_v2$PARAMCD]
  n_miss[subj_plot] <- sum(is.na(ymat))
  n_week[subj_plot] <- nrow(ymat)
}

treatment <- d_time[match(subj_unique, d_time$USUBJID), "newTRT"]
names(treatment) <- subj_unique

table(n_miss)
n_each_plot <- 4
n_week_keep <- 6
set.seed(1)
arms_plot <- c("Placebo", "300mg")
subj_plot_list <- list()
subj_plot_list[["Placebo"]] <- sample(subj_unique[treatment == "Placebo" & n_week == n_week_keep], size = n_each_plot)
subj_plot_list[["300mg"]] <- sample(subj_unique[treatment == "300mg" & n_week == n_week_keep], size = n_each_plot)

# subj_unique[which(n_miss == 8)[1]]#"CAIN457F2206_0153_00002"


# # control$col_arm["300mg"] <- control$col_arm["300mg"]
# graphics.off()
# export_plots <- TRUE
# if (export_plots){
#   pdf(file = paste0("figures/intro_to_factor_models.pdf"), 12, 8)
# }





# control$col_arm["300mg"] <- control$col_arm["300mg"]
graphics.off()
export_plots <- T
if (export_plots){
  pdf(file = paste0("figures/intro_to_stage_2.pdf"), 12, 5)
}



week_plot <- c(0, 2, 4, 8, 12, 16)
side_marg_wid <- .4
cexax <- .75
for (arm_plot in arms_plot[1:2]) {
  for (subj_num in 1:n_each_plot) {
    subj_plot <- subj_plot_list[[arm_plot]][subj_num]
    d_time_sub <- d_time[which(d_time$USUBJID == subj_plot & d_time$AVISITN %in% week_plot), ]
    week_plot_sub <- d_time_sub$AVISITN
    study <- d_time_sub$study[1]
    indication <- control$study_indic[grep(study, names(control$study_indic))]
    mfd_obj <- analysis_results_list[[paste0(indication, "_", study_plot)]]$mfd_object
    Z <- mfd_obj$Z[which(mfd_obj$data_out$USUBJID == subj_plot), ]
    S <- mfd_obj$S[which(mfd_obj$data_out$USUBJID == subj_plot), ]
    
    ymat <- d_time_sub[, nms_params_v2$PARAMCD]
    rownames(ymat) <- d_time_sub$AVISITN
    n_meas <- nrow(nms_params_v2)
    
    ##########################
    # Add scores
    par(mar = c(0, 0, 0, 0), oma = c(1, 1, 1, 1))
    fac_coords <- list()
    score_y_offset <- 0
    score_plot_hei <- .225
    mean_score_plot_wid <- .1
    mean_score_plot_hei <- .25
    
    centre_y <- .5
    treat_plac_offset <- .2
    score_plot_wid <- side_marg_wid / 1.5 / n_each_plot
    horiz_offset <- 0#score_plot_wid / 2
    for (fac_num in 1) {
      if (arm_plot == "Placebo") {
        fac_coords[[fac_num]] <- c(subj_num / (n_each_plot + 1) * side_marg_wid + c(-1, 1) * score_plot_wid / 2,
                                   centre_y + c(1, 0, -1)[fac_num] * score_y_offset + c(-1, 1) * score_plot_hei / 2 + treat_plac_offset)
      }
      if (arm_plot == "300mg") {
        fac_coords[[fac_num]] <- c(subj_num / (n_each_plot + 1) * side_marg_wid + c(-1, 1) * score_plot_wid / 2 + horiz_offset,
                                   centre_y + c(1, 0, -1)[fac_num] * score_y_offset + c(-1, 1) * score_plot_hei / 2 - treat_plac_offset)
      }
      suppressWarnings({
        par(fig = fac_coords[[fac_num]], new = T)
      })
      
      plot(week_plot, Z[, fac_num], ty = "n", 
           ylim = c(-1, 1) * 3,#base::range(c(c(Z + 2 * S), c(Z - 2 * S)), na.rm = T),
           xaxt = "n",
           yaxt = "n",
           las = 1,
           ylab = "",
           xlab = "")
      col_curr <- control$col_arm[arm_plot]
      lines(week_plot, Z[, fac_num], lwd = 2, col = col_curr)
      lines(week_plot, Z[, fac_num] + 2 * S[, fac_num], lwd = 1, col = col_curr)
      lines(week_plot, Z[, fac_num] - 2 * S[, fac_num], lwd = 1, col = col_curr)
      if (subj_num == 1) {
        axis(side = 2, las = 2, cex.axis = .8,
             at = c(-2, 0, 2), labels = c(-2, 0, 2))
      }
      axis(side = 1,
           at = c(0, 15),
           labels = NA,
           cex.axis = cexax,
           tcl = -.25)
      mtext(side = 1,
            at = c(0, 15),
            text = c(0, 15),
            cex = cexax)
      abline(h = 0, lty = 3)
      arm_lab <- ifelse(arm_plot == "Placebo", "Placebo", "300mg Treated")
      # lab_curr <- paste0(, " ", subj_num)
      subj_num_lab <- subj_num + n_each_plot * ifelse(arm_plot == "Placebo", 0, 1)
      mtext(side = 3, text = paste0("Patient ", subj_num_lab), cex = cexax, line = .5)
      mtext(side = 3, text = paste0("(", arm_lab, ")"), cex = cexax * .8, line = 0)
      # axis(side = 4, 
      #      cex.axis = cexax,
      #      las = 2)
      # mtext(side = 3, 
      #       text = paste0("Factor ", fac_num), 
      #       line = .15,
      #       cex = .9)
      # mtext(side = 1, 
      #       line = .05, 
      #       text = "Week", 
      #       cex = cexax)
    }
    # mtext(side = 3, text = "Outputs: scores, Z", line = 2)
  }  
}
# 
# for (fac_num in 1:3) {
#   mtext(outer = T,
#         at = c(5, 3, 1)[fac_num] / 6, 
#         side = 2, 
#         text = paste0("Factor ", fac_num),
#         las = 0)
# }

score_coords<-ri_coords<-ar1_coords<-res_coords <- list()
for (fac_num in 1) {
  score_coords[[fac_num]] <- c(side_marg_wid + (1 - side_marg_wid) * .175 + c(-1, 1) * mean_score_plot_wid / 2,
                               centre_y + c(1, 0, -1)[fac_num] * score_y_offset + c(-1, 1) * mean_score_plot_hei / 2)
  ri_coords[[fac_num]] <- c(side_marg_wid + (1 - side_marg_wid) * .45 + c(-1, 1) * mean_score_plot_wid / 2,
                            centre_y + c(1, 0, -1)[fac_num] * score_y_offset + c(-1, 1) * mean_score_plot_hei / 2)
  ar1_coords[[fac_num]] <- c(side_marg_wid + (1 - side_marg_wid) * .65 + c(-1, 1) * mean_score_plot_wid / 2,
                             centre_y + c(1, 0, -1)[fac_num] * score_y_offset + c(-1, 1) * mean_score_plot_hei / 2)
  res_coords[[fac_num]] <- c(side_marg_wid + (1 - side_marg_wid) * .85 + c(-1, 1) * mean_score_plot_wid / 2,
                             centre_y + c(1, 0, -1)[fac_num] * score_y_offset + c(-1, 1) * mean_score_plot_hei / 2)
  
  par(fig = score_coords[[fac_num]], new = T)
  plot_scores(control, 
              analysis_results_list,
              analysis_plot = analysis_to_plot,
              fac_plot = fac_num,
              estimate_type = "treatment_means",
              arms_in = arms_plot,
              add_legend = FALSE,
              cex.axis = .7,
              xaxt = "n",
              yaxt = "n",
              add_zero_horiz_line = F)
  # mtext(side = 3, text = paste0("Factor ", fac_num))
  abline(h = 0, lty = 3)
  axis(side = 2, las = 2, cex.axis = .8,
       at = c(-1, 0, 1), labels = c(-1, 0, 1))
  axis(side = 1,
       at = c(0, 15),
       labels = NA,
       cex.axis = cexax,
       tcl = -.25)
  mtext(side = 1,
        at = c(0, 15),
        text = c(0, 15),
        cex = cexax)
  tit_line1 <- 2
  tit_line2 <- 1.25
  mtext(side = 3, line = tit_line1, text = "Treatment")
  mtext(side = 3, line = tit_line2, text = "arm mean")
  estimate_type <- "estimate_type"
  legend_title <- switch(estimate_type,
                         estimate_type = "Means",
                         contrasts_vs_placebo = "Contrasts")
  # legend(x = max(t_unique) * 1.15, 
  #        # y = mean(par("usr")[3:4]), 
  #        y = par("usr")[3] + (par("usr")[4] - par("usr")[3]) * legend_y_mult,
   legend(x = "topright",
         legend = paste0(arms_plot, ifelse(estimate_type == "contrasts_vs_placebo", " - Placebo", "")), 
         col = 1,#control_plot$col_arm[arms_plot],
         pch = 22, 
         pt.bg = control$col_arm[arms_plot], 
         cex = .7,
         pt.cex = 1.5,
         yjust = .5
         # title = legend_title,
         # title.adj = 0.5
  )
  # mtext(side = 3, line = 4, text = "Linear predictor", cex = 1)

   
  outer_line <- 6 
  mtext(side = 3, line = outer_line, text = "(b) Gaussian mean", cex = 1)
  mtext(side = 3, line = outer_line - 1, text = "(Stage 2)", cex = 1)
  mtext(side = 3, line = outer_line, at = -40, text = "(a) Stage 1 output", cex = 1)
  # mtext(side = 3, line = outer_line, at = -5, text = "~ Normal(", cex = 1)
  ar1_covariance <- function(n_time, rho, noise_sd, intercept_sd) {
    R_AR1 <- rho ^ abs(outer(1:n_time, 1:n_time, "-")) # AR1 temporal correlation
    R_diagonal <- diag(n_time)
    R_ones <- matrix(1, n_time, n_time)
    delta_cov <- noise_sd^2 * R_AR1 + intercept_sd^2 * R_ones
    return(delta_cov)
  }
  
  
  ar1_cov <- ar1_covariance(n_time = 17, rho = .95, noise_sd = 1, intercept_sd = 0)
  ar1_cov <- ar1_cov[, rev(1:nrow(ar1_cov))]
  ri_cov <- ar1_cov
  ri_cov[] <- 1
  par(fig = ri_coords[[fac_num]], new = T)
  image(x = 1:17, y = 1:17, xaxt = "n", yaxt = "n", z = ri_cov, col = grey(0))
  axis(side = 1,
       at = c(1, 16),
       labels = NA,
       cex.axis = cexax,
       tcl = -.25)
  mtext(side = 1,
        at = c(1, 16),
        text = c(0, 15),
        cex = cexax)
  axis(side = 2,
       at = c(1, 16),
       labels = NA,
       cex.axis = cexax,
       tcl = -.25)
  mtext(side = 2,
        at = c(1, 16),
        text = c(0, 15),
        cex = cexax)
  mtext(side = 3, line = tit_line1, text = "Random")
  mtext(side = 3, line = tit_line2, text = "intercept")

  greys <- grey(seq(1, 0, len = 1000))
  par(fig = ar1_coords[[fac_num]], new = T)
  image(x = 1:17, y = 1:17, xaxt = "n", yaxt = "n", z = ar1_cov, col = greys)
  axis(side = 1,
       at = c(1, 16),
       labels = NA,
       cex.axis = cexax,
       tcl = -.25)
  mtext(side = 1,
        at = c(1, 16),
        text = c(0, 15),
        cex = cexax)
  axis(side = 2,
       at = c(1, 16),
       labels = NA,
       cex.axis = cexax,
       tcl = -.25)
  mtext(side = 2,
        at = c(1, 16),
        text = c(0, 15),
        cex = cexax)
  mtext(side = 3, line = tit_line1, text = "AR1")
  mtext(side = 3, line = tit_line2, text = "process")
  mtext(side = 3, line = mean(c(tit_line1, tit_line2)), at = -2, text = "+")
  
  
  mtext(side = 3, line = outer_line, text = "(c) Gaussian covariance structure                             ", cex = 1)
  mtext(side = 3, line = outer_line - 1, text = "(Stage 2)                             ", cex = 1)
  
  
  res_cov <- ar1_cov
  res_cov[] <- 0
  diag(res_cov) <- runif(17, .75, 1.25)
  res_cov <- res_cov[, rev(1:nrow(res_cov))]
  
  # par(fig = res_coords[[fac_num]], new = T)
  # image(x = 1:17, y = 1:17, xaxt = "n", yaxt = "n", z = res_cov, col = greys)
  # axis(side = 1,
  #      at = c(1, 16),
  #      labels = NA,
  #      cex.axis = cexax,
  #      tcl = -.25)
  # mtext(side = 1,
  #       at = c(1, 16),
  #       text = c(0, 15),
  #       cex = cexax)
  # axis(side = 2,
  #      at = c(1, 16),
  #      labels = NA,
  #      cex.axis = cexax,
  #      tcl = -.25)
  # mtext(side = 2,
  #       at = c(1, 16),
  #       text = c(0, 15),
  #       cex = cexax)
  # mtext(side = 3, line = tit_line1, text = "Diagonal")
  # mtext(side = 3, line = tit_line2, text = "residual")
  # mtext(side = 3, line = 0.5, text = "(Stage 1 outputted SDs)", cex = .75)
  # 
  # mtext(side = 3, line = mean(c(tit_line1, tit_line2)), at = -2, text = "+")
  # mtext(side = 3, line = 0, text = "covariance")
  
}
if (export_plots) {
  dev.off()
}




