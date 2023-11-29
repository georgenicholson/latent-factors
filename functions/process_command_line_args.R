process_command_line_args <- function(args_default = list(TASK_ID = 1)) {
  args_test_presence <- commandArgs(trailingOnly = TRUE)
  if (length(args_test_presence) > 0) {
    parser <- argparse::ArgumentParser()
    parser$add_argument("--TASK_ID", required = TRUE,
                        help = "Task number")
    args <- parser$parse_args()
    args$TASK_ID <- as.numeric(args$TASK_ID)
  } else {
    args <- args_default
  }
  return(args)
}
