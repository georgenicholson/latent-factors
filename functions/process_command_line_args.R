process_command_line_args <- function(args_default = list(TASK_ID = 1,
                                                          N_TASKS = 1000,
                                                          CONTROL_NUMBER = 0,
                                                          LINEAR_MODEL_FORM = "linear",
                                                          GEN_MODEL = "add")) {
  args_test_presence <- commandArgs(trailingOnly = TRUE)
  if (length(args_test_presence) > 0) {
    parser <- argparse::ArgumentParser()
    parser$add_argument("--n_subj", required = TRUE,
                        help = "Total number of subjects")
    parser$add_argument("--cor_ar1_lag", required = TRUE,
                        help = "Lag-1 autocorrelation in AR-1 process")
    parser$add_argument("--theta_RI_sd", required = TRUE,
                        help = "")
    parser$add_argument("--theta_resid_sd", required = TRUE,
                        help = "")
    parser$add_argument("--", required = TRUE,
                        help = "")
    parser$add_argument("--", required = TRUE,
                        help = "")
    parser$add_argument("--", required = TRUE,
                        help = "")
    parser$add_argument("--", required = TRUE,
                        help = "")
    parser$add_argument("--", required = TRUE,
                        help = "")
    parser$add_argument("--N_TASKS", required = TRUE,
                        help = "Number of tasks in array")
    parser$add_argument("--CONTROL_NUMBER", required = TRUE,
                        help = "Non-negative integer specifiying index of the analysis pipeline to run")
    parser$add_argument("--GEN_MODEL", required = FALSE,
                        help = "add/dom")
    parser$add_argument("--SEX_INTERACTION", required = FALSE,
                        help = "Genotype-sex interaction?")
    parser$add_argument("--LINEAR_MODEL_FORM", required = FALSE,
                        help = "constant/linear/quadratic")
    parser$add_argument("--CV_FOLD", required = FALSE,
                        help = "CV fold")
    parser$add_argument("--GEN_CONST", required = FALSE,
                        help = "Constant longitudinal genetic effects")
    parser$add_argument("--SNP_NUM", required = FALSE,
                        help = "Marker index")
    args <- parser$parse_args()
    args$TASK_ID <- as.numeric(args$TASK_ID)
    args$N_TASKS <- as.numeric(args$N_TASKS)
    args$CONTROL_NUMBER <- as.numeric(args$CONTROL_NUMBER)
    args$SEX_INTERACTION <- as.logical(as.numeric(args$SEX_INTERACTION))
    args$LINEAR_MODEL_FORM <- as.character(args$LINEAR_MODEL_FORM)
    args$CV_FOLD <- as.numeric(args$CV_FOLD)
    args$GEN_CONST <- as.logical(as.numeric(args$GEN_CONST))
    args$SNP_NUM <- as.numeric(args$SNP_NUM)
    
  } else {
    args <- args_default
  }
  return(args)
}




n_subj <- c(100)[1]
cor_ar1_lag <- c(.2, .75)[1]
theta_resid_sd<-theta_RI_sd <- c(.25, 1)[2]
endpoint_proportions <- list(c(0.6, 0.3, 0.1),
                             c(1, 1, 1, 1) / 4)[[2]]
missing_proportion <- c(0, .25)[1]
beta_Emax <- c(.1, .25, .5)[1]
y_resid_sd <- c(.25, 1, 2)[2]
