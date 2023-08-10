
# source("R/misc.R")
# source("R/mfd.R")

# data_to_conv$ord_cat <- factor(data_to_conv$ord_cat, levels = c("no response", "light response", "high response"), ordered=TRUE)

data_to_conv <- dplyr::as_tibble(data.frame(bin = factor(c("N", "Y", "N", "N", "Y"), levels=c("N", "Y")), 
                                            ord_cat = factor(c("no response", "high response", "light response", "light response", "no response"), levels = c("no response", "light response", "high response"), ordered=TRUE), 
                                            num=1:5, 
                                            char_col = c("a", "b", "c", "d", "e"), 
                                            unord_cat = factor(c("green", NA, "blue", "green", "red"), ordered=FALSE), 
                                            unord_cat_w_base = factor(c("no_event", "heartattack", "stroke", "incident", "stroke"), levels=c("no_event", "stroke", "incident", "heartattack"), ordered=FALSE )))

# cols <- c("bin", "ord_cat", "unord_cat", "unord_cat_w_base")
# col_types <- c("bin", "ord_cat", "unord_cat", "unord_cat_w_baseline")
# col_infos <- list(NULL, c("no response", "light response", "high response"), NA, "no_event")

# convert_factor_to_numeric(data=data_to_conv)

test_that("'functionality': checking all possible col_types once", {
  expect_error(convert_factor_to_numeric(data=data_to_conv), NA)
})

test_that("input: checking all possible col_types once", {
  expect_error(convert_factor_to_numeric(data=as.data.frame(data_to_conv)), "data is not a tibble.")
})

