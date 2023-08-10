

sample_units <- c()
for (i in 1:15) {
  sample_units <- c(sample_units, rep(paste0("subj", as.character(i)), 2))
}
data_test <- tibble(a = 1:30, indic = c(rep("a", 10), rep("b", 20)), subj_ID = sample_units )

test_that("no error", {
  expect_error(sampling(data = data_test, series_group = "subj_ID", subgroup = "indic", total_sample_size = 5, type = "equal_n_series_from_subgroups"), NA)
})

test_that("input: data is a tibble", {
  expect_error(sampling(data = as.data.frame(data_test), series_group = "subj_ID", subgroup = "indic", total_sample_size = 5, type = "equal_n_series_from_subgroups"), "data is not a tibble.")
})

test_that("input: total_sample_size is a positive integer", {
  expect_error(sampling(data = data_test, series_group = "subj_ID", subgroup = "indic", total_sample_size = 5.5, type = "equal_n_series_from_subgroups"), "total_sample_size must be positive int.")
  expect_error(sampling(data = data_test, series_group = "subj_ID", subgroup = "indic", total_sample_size = -2, type = "equal_n_series_from_subgroups"), "total_sample_size must be positive int.")
})

test_that("input: type is valid", {
  expect_error(sampling(data = data_test, series_group = "subj_ID", subgroup = "indic", total_sample_size = 5, type = "something else"), "type specification is not valid.")
})

test_that("input: series_group and subgroup columns must be present in data", {
  expect_error(sampling(data = data_test %>% select(-subj_ID), series_group = "subj_ID", subgroup = "indic", total_sample_size = 5, type = "equal_n_series_from_subgroups"), "Columns series_group or subgroup not present in data.")
  expect_error(sampling(data = data_test %>% select(-indic), series_group = "subj_ID", subgroup = "indic", total_sample_size = 5, type = "equal_n_series_from_subgroups"), "Columns series_group or subgroup not present in data.")
})

# TODO functionality tests
