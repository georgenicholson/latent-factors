source("scripts/04_plotting_intro.R")
########################################
# Plot single study results
include_loadings <- FALSE
# export_plots <- T#FALSE
for (indication in c("PsA", "RA")) {
  for (estimate_type in score_plot_types) {
    if (export_plots){
      pdf(file = paste0("figures/meta_analysis_", estimate_type, "_indication_", indication, ".pdf"), 12, 10)
    }
    plot_num <- 0
    letter_cex <- 1.25
    if (include_loadings) {
      layout(mat = matrix(c(1, 2, 3, 4, 1, 5, 6, 7, 1, 8, 9, 10), nrow = 3, ncol = control$n_fac + 1, byrow = TRUE))
      par(mar = c(2, 2, 2, 2), oma = c(4, 10, 5, 17))
      plot_loadings(control, 
                    analysis_results_list, 
                    analysis_plot = "PsA_F2312_cross_indication_loadings",
                    add_legend = T)
      add_letter(letters[plot_num <- plot_num + 1], y_fac = .025, cex = letter_cex)
    } else {
      par(mar = c(2, 2, 2, 2), 
          oma = c(4, 4, 5, 17),
          mfrow = c(3, 3))
    }
    # indication_meta <- list(PsA = list(studies_in = c("PsA_F2312_cross_indication_loadings",
    #                                 "PsA_F2342_cross_indication_loadings",
    #                                 "PsA_meta_cross_indication_loadings"),
    #                                 arms_in = c("300mg", "150mg")),
    #                         RA = list(studies_in = c("RA_F2302_cross_indication_loadings",
    #                                                  "RA_F2309_cross_indication_loadings",
    #                                                  "RA_meta_cross_indication_loadings"),
    #                                   arms_in = c("75mg", "150mg")))
    indication_meta <- list(PsA = list(studies_in = c("PsA_F2312",
                                                      "PsA_F2342",
                                                      "PsA_meta"),
                                       arms_in = c("300mg", "150mg")),
                            RA = list(studies_in = c("RA_F2302",
                                                     "RA_F2309",
                                                     "RA_meta"),
                                      arms_in = c("75mg", "150mg")))
    studies_to_plot <- indication_meta[[indication]]$studies_in
    for (study_to_plot in studies_to_plot) {
      study_lab <- paste0(indication, ", ", control$do_list[[study_to_plot]]$stage2, "")
      if (grepl(indication, study_to_plot)) {
        study_lab <- paste0(indication, ", ", paste(control$do_list[[study_to_plot]]$stage2, collapse = " & "), "")
      }
      for (fac_curr in 1:control$n_fac) {
        plot_scores(control, 
                    analysis_results_list,
                    analysis_plot = study_to_plot,
                    fac_plot = fac_curr,
                    estimate_type = estimate_type,
                    add_legend = fac_curr == control$n_fac & match(study_to_plot, indication_meta[[indication]]$studies_in) == 1,
                    arms_in = indication_meta[[indication]]$arms_in)
        plotnum <- plot_num + 1
        if (fac_curr == 1) {
          add_letter(letters[plot_num <- plot_num + 1], 
                     y_fac = .1,
                     x_fac = .1)
        }
        if (study_to_plot == studies_to_plot[1]) {
          mtext(side = 3, text = control$fac_labs[fac_curr], line = 2)
        }
      }
      if (fac_curr == control$n_fac) {
        add_study_lab(control = control, 
                      study_to_plot = study_to_plot)
        # mtext(side = 4, text = study_lab, cex = 1.2, las = 2, line = 2)
      }
    }
    mtext(outer = TRUE,
          side = 1,
          text = "Week",
          line = 1)
    if (export_plots){
      dev.off()
    }
  }
}
  
#   layout(mat = matrix(c(1, 1, 1, 2, 3, 4), nrow = 3, ncol = 3, byrow = FALSE))
#   par(mfcol = c(4, 4), mar = c(3, 3, 3, 3), oma = c(4, 10, 5, 10))
#   
#   single_studies <- c("F2312", "F2342", "F2302", "F2309")
#   for (study_to_plot in single_studies) {
#     plot_loadings(control, 
#                   analysis_results_list, 
#                   analysis_plot = study_to_plot,
#                   add_legend = study_to_plot == single_studies[1])
#     add_letter(letters[plot_num <- plot_num + 1], y_fac = .025)
#   }
#   for (fac_curr in 1:control$n_fac) {
#     for (study_to_plot in single_studies) {
#       plot_scores(control, 
#                   analysis_results_list,
#                   analysis_plot = study_to_plot,
#                   fac_plot = fac_curr,
#                   estimate_type = estimate_type,
#                   add_legend = fac_curr == control$n_fac)
#       plotnum <- plot_num + 1
#       add_letter(letters[plot_num <- plot_num + 1])
#     }
#   }
#   mtext(outer = TRUE, at = seq(.25, .75, by = .25) + .125, text = control$fac_labs)
#   if (export_plots){
#     dev.off()
#   }
#   
# }

