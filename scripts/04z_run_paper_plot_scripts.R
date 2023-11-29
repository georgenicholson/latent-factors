
#######################################
# Generate plots
export_plots <- TRUE
source("scripts/04b_plot_single_study.R")
source("scripts/04d_meta_analysis.R")
source("scripts/04f_compare_indications.R")
source("scripts/04f_compare_indications_same_plot.R")
source("scripts/04g_compare_loadings.R")
source("scripts/04h_intro_to_factor_analysis_plot.R")
source("scripts/04i_intro_to_stage_2_single_factor.R")
system("rm plots.tar.gz")
system("(cd figures && tar cvzf ../plots.tar.gz *)")
browseURL("plots.tar.gz")


dir.create("output/results_out", showWarnings = FALSE)
file.copy(from = "output/analysis_results_list.RDS", to = "output/results_out/analysis_results_list.RDS")
system("(cd output/results_out && tar cvzf ../results_out.tar.gz *)")
browseURL("output/results_out.tar.gz")


