#######################################
# Source functions
fn_files_to_source <- list.files("functions/", full.names = TRUE)
for (file_curr in fn_files_to_source) {
  source(file_curr)
}

control <- get_control_parameters_simulations(parameter_set = 1)
######################################################
# Export the full set of tables and emax plot in tar.gz format
source("scripts/s01_emax_curves_demo_plot.R")
source("scripts/s03_output_latex_tables.R")
source("scripts/s05_output_benchmarking_latex_tables.R")
tab_export_files <- c(control$default_sim_table_file, 
                      control$sim_results_file, 
                      control$emax_example_plot_file, 
                      control$numbers_out_dir,
                      control$benchmarking_results_files,
                      control$benchmarking_runtime_file)
tab_export_files_zipped <- zip::zip(zipfile = file.path("zipped_tables.zip"),
                                    files = tab_export_files,
                                    mode = "cherry-pick")
tab_export_files_zipped <- file.path("zipped_tables.tar.gz")
tar(tarfile = tab_export_files_zipped,
    files = tab_export_files,
    compression = 'gzip',
    tar = "tar")
browseURL(tab_export_files_zipped) # this downloads locally to Downloads folder on a PC


