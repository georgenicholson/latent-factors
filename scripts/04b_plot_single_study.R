source("scripts/04_plotting_intro.R")

########################################
# Plot single study results
study_to_plot <- "F2309"
for (study_to_plot in control$single_studies) {
  if (export_plots){
    pdf(file = paste0("figures/single_study_", study_to_plot, ".pdf"), 12, 7)
  }
  plot_num <- 0
  layout(mat = matrix(c(1, 2, 3, 4, 1, 5, 6, 7), nrow = 2, ncol = control$n_fac + 1, byrow = TRUE))
  par(mar = c(3, 3, 3, 3), oma = c(4, 10, 5, 15))
  analysis_plot <- study_to_plot
  plot_loadings(control, 
                analysis_results_list, 
                analysis_plot = study_to_plot,
                add_legend = T)
  add_letter(letters[plot_num <- plot_num + 1], y_fac = .025)
  par(new = F)
  for (estimate_type in score_plot_types) {
    for (fac_curr in 1:control$n_fac) {
        plot_scores(control, 
                  analysis_results_list,
                  analysis_plot = study_to_plot,
                  fac_plot = fac_curr,
                  estimate_type = estimate_type,
                  add_legend = fac_curr == control$n_fac)
      plotnum <- plot_num + 1
      if (fac_curr == 1) {
        add_letter(letters[plot_num <- plot_num + 1])
      }
      if (fac_curr == control$n_fac) {
        add_study_lab(control = control, 
                      study_to_plot = study_to_plot)
      }
    }
  }
  mtext(outer = TRUE, at = seq(.25, .75, by = .25) + .125, text = control$fac_labs)
  if (export_plots){
    dev.off()
  }
}

