process_command_line_args <- function(args_default = list(TASK_ID = 1,
                                                          N_TASKS = 1000)) {
  args_test_presence <- commandArgs(trailingOnly = TRUE)
  if (length(args_test_presence) > 0) {
    parser <- argparse::ArgumentParser()
    parser$add_argument("--TASK_ID", required = TRUE,
                        help = "Task number")
    parser$add_argument("--N_TASKS", required = TRUE,
                        help = "Number of tasks in array")
    args <- parser$parse_args()
    args$TASK_ID <- as.numeric(args$TASK_ID)
    args$N_TASKS <- as.numeric(args$N_TASKS)
  } else {
    args <- args_default
  }
  return(args)
}
