test_that("input: data is tibble", {
  sample_units <- c()
  for (i in 1:15) {
    sample_units <- c(sample_units, rep(paste0("subj", as.character(i)), 2))
  }
  expect_error(stratified_sampling(data = data.frame(a = 1:30, b = c(rep("a", 10), rep("b", 20)),
                                                    c = sample_units ),
                                   group = "b", total_sample_size = 10,
                                   sample_unit = "c", with_replace = FALSE)  )
})

test_that("input: group variable in data", {
  sample_units <- c()
  for (i in 1:15) {
    sample_units <- c(sample_units, rep(paste0("subj", as.character(i)), 2))
  }
  expect_error(stratified_sampling(data = tibble(a = 1:30, b = c(rep("a", 10), rep("b", 20)),
                                                        c = sample_units ),
                                   group = "d", total_sample_size = 10,
                                   sample_unit = "c", with_replace = FALSE)  )
})

test_that("input: total_sample_size is positive integer", {
  sample_units <- c()
  for (i in 1:15) {
    sample_units <- c(sample_units, rep(paste0("subj", as.character(i)), 2))
  }
  expect_error(stratified_sampling(data = tibble(a = 1:30, b = c(rep("a", 10), rep("b", 20)),
                                                    c = sample_units ),
                                   group = "b", total_sample_size = 10.1,
                                   sample_unit = "c", with_replace = FALSE)  )
  expect_error(stratified_sampling(data = tibble(a = 1:30, b = c(rep("a", 10), rep("b", 20)),
                                                    c = sample_units ),
                                   group = "b", total_sample_size = -3,
                                   sample_unit = "c", with_replace = FALSE)  )
})


test_that("input: sample_unit variable in columns of data", {
  sample_units <- c()
  for (i in 1:15) {
    sample_units <- c(sample_units, rep(paste0("subj", as.character(i)), 2))
  }
  expect_error(stratified_sampling(data = tibble(a = 1:30, b = c(rep("a", 10), rep("b", 20)),
                                                    c = sample_units ),
                                   group = "b", total_sample_size = 10,
                                   sample_unit = "d", with_replace = FALSE)  )
})


test_that("input: subgroup_sample_strat correctly specified", {
  sample_units <- c()
  for (i in 1:15) {
    sample_units <- c(sample_units, rep(paste0("subj", as.character(i)), 2))
  }
  expect_error(stratified_sampling(data = tibble(a = 1:30, b = c(rep("a", 10), rep("b", 20)),
                                                    c = sample_units ),
                                   group = "b", total_sample_size = 10,
                                   sample_unit = "c", subgroup_sample_strat = "something else", with_replace = FALSE)  )
})


test_that("input: with_replace is boolean", {
  sample_units <- c()
  for (i in 1:15) {
    sample_units <- c(sample_units, rep(paste0("subj", as.character(i)), 2))
  }
  expect_error(stratified_sampling(data = tibble(a = 1:30, b = c(rep("a", 10), rep("b", 20)),
                                                    c = sample_units ),
                                   group = "b", total_sample_size = 10,
                                   sample_unit = "c", with_replace = "no boolean")  )
})


test_that("input: total_sample_size too large", {
  sample_units <- c()
  for (i in 1:15) {
    sample_units <- c(sample_units, rep(paste0("subj", as.character(i)), 2))
  }
  expect_error(stratified_sampling(data = tibble(a = 1:30, b = c(rep("a", 10), rep("b", 20)),
                                                    c = sample_units ),
                                   group = "b", total_sample_size = 16,
                                   sample_unit = "c", with_replace = FALSE)  )
})


test_that("functionality: sampling ratio roughly ok", {
  sample_units <- c()
  for (i in 1:15000) {
    sample_units <- c(sample_units, rep(paste0("subj", as.character(i)), 2))
  }
  data <- tibble(a = 1:30000, b = c(rep("a", 10000), rep("b", 20000)),
                    c = sample_units )
  data_sampled <- stratified_sampling(data = data,
                                   group = "b", total_sample_size = 10000,
                                   sample_unit = "c", with_replace = FALSE)
  a_samples <- data %>% filter(b == "a") %>% select(c) %>% n_distinct()
  b_samples <- data %>% filter(b == "b") %>% select(c) %>% n_distinct()
  expect_equal(round(a_samples / (a_samples + b_samples), digits=2), 0.33)
  expect_equal(round(b_samples / (a_samples + b_samples), digits=2), 0.67)
})


# TODO test grouping works correctly




