source("scripts/04_plotting_intro.R")
########################################
# Plot single study results
# export_plots <- F
if (export_plots){
  pdf(file = paste0("figures/compare_loadings.pdf"), 12, 8)
}
plot_num <- 0

plot_pooled_loadings <- FALSE
layout(mat = matrix(c(1, 2, 3, 4, 1, 5, 6, 7, 1, 8, 9, 10, 1, 11, 12, 13), nrow = 4, ncol = control$n_fac + 1, byrow = TRUE))
par(mar = c(2, 1, 1, 2), 
    oma = c(5, 8, 7, 2),
    mfrow = c(1, 6 + ifelse (plot_pooled_loadings, 1, 0)))

names(control$do_list)
analyses_to_plot <- c("PsA_F2312", "PsA_F2342", "PsA_meta", "RA_F2302", "RA_F2309", "RA_meta")
if (plot_pooled_loadings) {
  analyses_to_plot <- c(analyses_to_plot,  "RA_F2302_cross_indication_loadings")
}
for (analysis_plot in analyses_to_plot) {
  plot_loadings(control, 
                analysis_results_list, 
                analysis_plot = analysis_plot,
                add_legend = analysis_plot == analyses_to_plot[3],
                add_endpt_labels = analysis_plot == analyses_to_plot[1],
                y_at_fac = .09)
  # plotnum <- plot_num + 1
  # add_letter(letters[plot_num <- plot_num + 1], 
  #            y_fac = 0.034,
  #            cex = 1.25)
  if (analysis_plot == "RA_F2302_cross_indication_loadings") {
    mtext(side = 3, text = "Pooled", cex = 1, line = 1) 
  } else {
    add_study_lab_top(control = control, 
                      study_to_plot = analysis_plot)
  }
  mtext(side = 3,
        at = 0.45,
        text = paste0("(", letters[match(analysis_plot, analyses_to_plot)], ")"),
        line = 1.75)
}


# 
# 
# add_letter(letters[plot_num <- plot_num + 1], y_fac = .025)
# par(new = F)
# # for (estimate_type in score_plot_types) {
# studies_to_plot <- c("PsA_meta_cross_indication_loadings", 
#                      "RA_meta_cross_indication_loadings")
# for (estimate_type in score_plot_types) {
#   if (estimate_type == "treatment_means") {
#     arms_in_curr <- c("Placebo" , "150mg")
#   }
#   if (estimate_type == "contrasts_vs_placebo") {
#     arms_in_curr <- c("150mg")
#   }
#   for (study_to_plot in studies_to_plot) {
#     for (fac_curr in 1:control$n_fac) {
#       plot_scores(control, 
#                   analysis_results_list,
#                   analysis_plot = study_to_plot,
#                   fac_plot = fac_curr,
#                   estimate_type = estimate_type,
#                   # add_legend = fac_curr == control$n_fac & match(study_to_plot, studies_to_plot) == 2,
#                   add_legend = fac_curr == control$n_fac & study_to_plot == studies_to_plot[2],
#                   arms_in = arms_in_curr,
#                   legend_y_mult = 1.225)
#       plotnum <- plot_num + 1
#       add_letter(letters[plot_num <- plot_num + 1], y_fac = 0.1)
#       if (fac_curr == control$n_fac) {
#         add_study_lab(study_to_plot = study_to_plot, eps_lab = .05)
#       }
#     }
#   }
# }
# mtext(outer = TRUE, at = seq(.25, .75, by = .25) + .125, text = control$fac_labs)
if (export_plots){
  dev.off()
}
