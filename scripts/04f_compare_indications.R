source("scripts/04_plotting_intro.R")
########################################
# Plot single study results
include_loadings <- FALSE
# export_plots <- FALSE
if (export_plots){
  pdf(file = paste0("figures/cross_indication_means_and_contrasts.pdf"), 12, 12)
}
plot_num <- 0
letter_cex <- 1.25
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
  par(mar = c(2, 2, 2, 2), 
      oma = c(4, 4, 5, 17),
      mfrow = c(4, 3))
}
par(new = F)
# for (estimate_type in score_plot_types) {
# studies_to_plot <- c("PsA_meta_cross_indication_loadings", 
#                      "RA_meta_cross_indication_loadings")
studies_to_plot <- c("PsA_meta", 
                     "RA_meta")
for (estimate_type in score_plot_types) {
  if (estimate_type == "treatment_means") {
    arms_in_curr <- c("Placebo" , "150mg")
  }
  if (estimate_type == "contrasts_vs_placebo") {
    arms_in_curr <- c("150mg")
  }
  for (study_to_plot in studies_to_plot) {
    for (fac_curr in 1:control$n_fac) {
      plot_scores(control, 
                  analysis_results_list,
                  analysis_plot = study_to_plot,
                  fac_plot = fac_curr,
                  estimate_type = estimate_type,
                  # add_legend = fac_curr == control$n_fac & match(study_to_plot, studies_to_plot) == 2,
                  add_legend = fac_curr == control$n_fac & study_to_plot == studies_to_plot[2],
                  arms_in = arms_in_curr,
                  legend_y_mult = 1.225)
      plotnum <- plot_num + 1
      add_letter(letters[plot_num <- plot_num + 1], y_fac = 0.1)
      if (fac_curr == control$n_fac) {
        add_study_lab(control = control, 
                      study_to_plot = study_to_plot, eps_lab = .075)
      }
      if (study_to_plot == studies_to_plot[1] & estimate_type == score_plot_types[1]) {
        mtext(side = 3, text = control$fac_labs[fac_curr], line = 2)
      }
    }
  }
}
mtext(outer = TRUE,
      side = 1,
      text = "Week",
      line = 1)
# mtext(outer = TRUE, at = seq(.25, .75, by = .25) + .125, text = control$fac_labs)
if (export_plots){
  dev.off()
}

# 
# 
# #######################################################################
# # Older version with separate plots for means and contrasts
# for (estimate_type in score_plot_types) {
#   if (estimate_type == "treatment_means") {
#     arms_in_curr <- c("Placebo" , "150mg")
#   }
#   if (estimate_type == "contrasts_vs_placebo") {
#     arms_in_curr <- c("150mg")
#   }
#   export_plots <- T
#   if (export_plots){
#     pdf(file = paste0("figures/cross_indication_",  estimate_type, ".pdf"), 12, 8)
#   }
#   plot_num <- 0
#   plot_num <- 0
#   layout(mat = matrix(c(1, 2, 3, 4, 1, 5, 6, 7), nrow = 2, ncol = control$n_fac + 1, byrow = TRUE))
#   par(mar = c(2, 2, 2, 2), oma = c(4, 10, 5, 15))
#   plot_loadings(control, 
#                 analysis_results_list, 
#                 analysis_plot = "F2312_cross_indication_loadings",
#                 add_legend = T)
#   add_letter(letters[plot_num <- plot_num + 1], y_fac = .025)
#   par(new = F)
#   # for (estimate_type in score_plot_types) {
#   studies_to_plot <- c("PsA_meta_cross_indication_loadings", 
#                        "RA_meta_cross_indication_loadings")
#   for (study_to_plot in studies_to_plot) {
#     indication <- ifelse(grepl("RA", study_to_plot), "RA", "PsA")
#     for (fac_curr in 1:control$n_fac) {
#       plot_scores(control, 
#                   analysis_results_list,
#                   analysis_plot = study_to_plot,
#                   fac_plot = fac_curr,
#                   estimate_type = estimate_type,
#                   add_legend = fac_curr == control$n_fac & match(study_to_plot, studies_to_plot) == 1,
#                   arms_in = arms_in_curr)
#       plotnum <- plot_num + 1
#       add_letter(letters[plot_num <- plot_num + 1])
#       if (fac_curr == control$n_fac) {
#         mtext(side = 4, text = indication, cex = 1.2, las = 2, line = 2)
#       }
#     }
#   }
#   mtext(outer = TRUE, at = seq(.25, .75, by = .25) + .125, text = control$fac_labs)
#   if (export_plots){
#     dev.off()
#   }
# }
