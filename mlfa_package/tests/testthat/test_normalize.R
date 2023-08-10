
data_test <- dplyr::tibble(a = 1:3, b = 4:6, c = 7:9, subj_ID = c("a", "a", "b"))

test_that("no error", {
  expect_error(normalize(data = data_test, columns = c("a", "b"), type = "per_group_center_per_meas_standard", series_group = "subj_ID"), NA)
  expect_error(normalize(data = data_test, columns = c("a", "b"), type = "per_meas_standard", series_group = "subj_ID"), NA)
})

test_that("input: data must be a tibble", {
  expect_error(normalize(data = as.data.frame(data_test), columns = c("a", "b"), type = "per_group_center_per_meas_standard", series_group = "subj_ID"), "data is not a tibble.")
})

test_that("input: columns must be in data", {
  expect_error(normalize(data = data_test, columns = c("a", "d"), type = "per_group_center_per_meas_standard", series_group = "subj_ID"), "Variable columns not in columns of data.")
})

test_that("type: must be one of the ones specified", {
  expect_error(normalize(data = data_test, columns = c("a", "b"), type = "something else", series_group = "subj_ID"), "Variable type incorrectly specified.")
})

test_that("input: column series_group must be present in data in type per_group_center_per_meas_standard", {
  expect_error(normalize(data = select(data_test, -subj_ID), columns = c("a", "b"), type = "per_group_center_per_meas_standard", series_group = "something else"), "series_group column not present in data.")
  # if per_meas_standard: series_group column need not be specified
  expect_error(normalize(data = select(data_test, -subj_ID), columns = c("a", "b"), type = "per_meas_standard", series_group = NULL), NA)
})


# TODO functionality tests
