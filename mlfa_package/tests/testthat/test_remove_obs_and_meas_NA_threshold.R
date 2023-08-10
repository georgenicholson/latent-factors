test_that("input: data must be tibble", {
  expect_error(remove_obs_and_meas_NA_threshold(data = data.frame(a = c(1, NA, 3), b = c(1, NA, 3), c = c(1, 2, 3)), columns = c("c"),
                                                min_nonNA_fraction_per_row = 1.0, min_nonNA_fraction_per_col = 0.0), "data is not a tibble.")
})

test_that("input: columns must be available AND not NULL", {
  expect_error(remove_obs_and_meas_NA_threshold(data = tibble(a = c(1, NA, 3), b = c(1, NA, 3), c = c(1, 2, 3)), columns = c("d"),
                                   min_nonNA_fraction_per_row = 1.0, min_nonNA_fraction_per_col = 0.0),
               "Variable columns specifies columns in data that are not existent.")
  expect_error(remove_obs_and_meas_NA_threshold(data = tibble(a = c(1, NA, 3), b = c(1, NA, 3), c = c(1, 2, 3)), columns = NULL,
                                                min_nonNA_fraction_per_row = 1.0, min_nonNA_fraction_per_col = 0.0),
               "Variable columns must not be NULL.")
})

test_that("input: columns must be available", {
  expect_error(remove_obs_and_meas_NA_threshold(data = tibble(a = c(1, NA, 3), b = c(1, NA, 3), c = c(1, 2, 3)), columns = c("a"),
                                                min_nonNA_fraction_per_row = 2.0, min_nonNA_fraction_per_col = 0.0))
  expect_error(remove_obs_and_meas_NA_threshold(data = tibble(a = c(1, NA, 3), b = c(1, NA, 3), c = c(1, 2, 3)), columns = c("a"),
                                                min_nonNA_fraction_per_row = -0.5, min_nonNA_fraction_per_col = 0.0))
  expect_error(remove_obs_and_meas_NA_threshold(data = tibble(a = c(1, NA, 3), b = c(1, NA, 3), c = c(1, 2, 3)), columns = c("a"),
                                                min_nonNA_fraction_per_row = 1.0, min_nonNA_fraction_per_col = 1.5))
  expect_error(remove_obs_and_meas_NA_threshold(data = tibble(a = c(1, NA, 3), b = c(1, NA, 3), c = c(1, 2, 3)), columns = c("a"),
                                                min_nonNA_fraction_per_row = 1.0, min_nonNA_fraction_per_col = -0.5))
})

test_that("functionality: remove_obs_and_meas_NA_threshold correctly removes ROWS", {
  data <- tibble(a = 1:3, b = c(1, NA, 3))
  expect_equal(remove_obs_and_meas_NA_threshold(data = tibble(a = c(1, NA, 3), b = c(1, NA, 3), c = c(1, 2, 3)), columns = c("a", "b", "c"),
                                                min_nonNA_fraction_per_row = 1.0, min_nonNA_fraction_per_col = 0.0),
                                                tibble(a = c(1, 3), b = c(1, 3), c = c(1, 3)))
})

test_that("functionality: remove_obs_and_meas_NA_threshold correctly removes COLUMNS", {
  data <- tibble(a = 1:3, b = c(1, NA, 3))
  expect_equal(remove_obs_and_meas_NA_threshold(data = tibble(a = c(1, NA, 3), b = c(1, NA, 3), c = c(1, 2, 3)), columns = c("a", "b", "c"),
                                                min_nonNA_fraction_per_row = 0.0, min_nonNA_fraction_per_col = 1.0), tibble(c = c(1, 2, 3)))
})

test_that("functionality: remove_obs_and_meas_NA_threshold correctly removes ROWS, only considering selected columns", {
  data <- tibble(a = 1:3, b = c(1, NA, 3))
  expect_equal(remove_obs_and_meas_NA_threshold(data = tibble(a = c(1, NA, 3), b = c(1, NA, 3), c = c(1, 2, 3)), columns = c("c"),
                                                               min_nonNA_fraction_per_row = 1.0, min_nonNA_fraction_per_col = 1.0),
                              tibble(a = c(1, NA, 3), b = c(1, NA, 3), c = c(1, 2, 3)) )  # do not remove anything
})


