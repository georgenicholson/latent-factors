#############################################
# ------------------------------------------
# Set up directory structure of repository:
# ------------------------------------------
#############################################

get_DIR <- function() {
  DIR <- list(data = "/data/as/processed/clinical/database_views/clinical_to_replace",
            out = "output",
            out_tmp_DATA = file.path("/data/users", Sys.info()["user"], "tmp_DATA_project1"),
            out_RA = "output/out_RA",
            out_PsA = "output/out_PsA",
            out_PsA_induction = "output/out_PsA/Induction",
            out_PsA_year1 = "output/out_PsA/Year1",
            out_pooled_RA_PsA = "output/out_pooled_RAPsA",
            meta = "meta",
            path_package = "mlfa_package/")
  
  DIR
}


