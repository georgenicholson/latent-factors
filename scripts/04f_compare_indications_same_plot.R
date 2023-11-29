source("scripts/04_plotting_intro.R")
include_loadings <- FALSE
# export_plots <- FALSE
if (export_plots){
  pdf(file = paste0("figures/cross_indication_means_and_contrasts_same_plot.pdf"), 12, 8)
}
plot_num <- 0
letter_cex <- 1.25
t_unique <- as.numeric(rownames(analysis_results_list[[1]][[1]][[1]]$m))
if (include_loadings) {
  layout(mat = matrix(c(1, 2, 3, 4, 1, 5, 6, 7, 1, 8, 9, 10, 1, 11, 12, 13), 
                      nrow = 4, 
                      ncol = control$n_fac + 1, 
                      byrow = TRUE)
  )
  par(mar = c(2, 2, 2, 2), 
      oma = c(4, 7, 5, 16))
  plot_loadings(control, 
                analysis_results_list, 
                analysis_plot = "PsA_F2312_cross_indication_loadings",
                add_legend = T)
  add_letter(letters[plot_num <- plot_num + 1], 
             y_fac = .025, 
             cex = letter_cex)
  layout(mat = matrix(c(1, 2, 3, 4, 1, 5, 6, 7, 1, 8, 9, 10), nrow = 3, ncol = control$n_fac + 1, byrow = TRUE))
} else {
  par(mar = c(2, 2, 4, 2), 
      oma = c(4, 4.5, 2, 17),
      mfrow = c(2, 3))
}
par(new = F)
studies_to_plot <- c("PsA_meta_cross_indication_loadings", 
                     "RA_meta_cross_indication_loadings")
for (estimate_type in score_plot_types) {
  if (estimate_type == "treatment_means") {
    arms_in_curr <- c("Placebo" , "150mg")
  }
  if (estimate_type == "contrasts_vs_placebo") {
    arms_in_curr <- c("150mg")
  }
  for (fac_curr in 1:control$n_fac) {
    plot_scores(control, 
                analysis_results_list,
                analysis_plot = studies_to_plot,
                fac_plot = fac_curr,
                estimate_type = estimate_type,
                add_legend = fac_curr == control$n_fac,
                arms_in = arms_in_curr,
                legend_y_mult = .5)
    # plotnum <- plot_num + 1
    if (fac_curr == 1) {
      add_letter(letters[plot_num <- plot_num + 1], 
                 y_fac = 0.08, 
                 x_fac = 0.09,
                 cex = 1.4)
    }
    if (fac_curr == 3 & estimate_type == score_plot_types[1]) {
      indications <- ifelse(grepl("RA", studies_to_plot), "RA", "PsA")
      studies <- sapply(control$do_list[studies_to_plot], function(x) paste0(paste(x$stage2, collapse = " & ")))
      study_lab <- indications#paste0(indications, " (", studies, ")")
      
      par(xpd = NA)
      legend(x = max(t_unique) * 1.15, 
             y = par("usr")[3] + (par("usr")[4] - par("usr")[3]) * 1.11,
             legend = study_lab, 
             lty = 1:length(studies_to_plot),
             lwd = 2,
             cex = 1.5,
             yjust = .5,
             # xjust = .5,
             title = "Indications",
             title.adj = 0.5
      )
      
      par(xpd = FALSE)
    }
    if (estimate_type == score_plot_types[1]) {
      mtext(side = 3, text = control$fac_labs[fac_curr], line = 2)
    }
  }
}
mtext(outer = TRUE,
      side = 1,
      text = "Week",
      line = 1)
mtext(outer = TRUE,
      side = 2,
      text = "Factor score",
      line = 2.5)
if (export_plots){
  dev.off()
}

