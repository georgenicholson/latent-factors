
Y_loadings <- matrix(1, 5, 5)
Y_scores <- matrix(3, 5, 5)

data_loadings <- cbind(Y_loadings, matrix(2, 5, 1))
colnames(data_loadings) <- c(1:6)
data_loadings <- as_tibble(data_loadings)

data_scores <- cbind(Y_scores, matrix(4, 5, 1))
colnames(data_scores) <- c(1:6)
data_scores <- as_tibble(data_scores)

test_that("input: no error", {
  expect_error(remove_all_NA_rows(Y_loadings, Y_scores, data_loadings, data_scores), NA)
})

test_that("input: Y_loadings is not a matrix", {
  expect_error(remove_all_NA_rows(as_tibble(Y_loadings), Y_scores, data_loadings, data_scores), "Y_loadings is not a matrix.")
})

test_that("input: Y_sores is not a matrix", {
  expect_error(remove_all_NA_rows(Y_loadings, as_tibble(Y_scores), data_loadings, data_scores), "Y_scores is not a matrix.")
})

test_that("input: data_loadings is not a tibble", {
  expect_error(remove_all_NA_rows(Y_loadings, Y_scores, as.matrix(data_loadings), data_scores), "data_loadings is not a tibble")
})

test_that("input: data_scores is not a tibble", {
  expect_error(remove_all_NA_rows(Y_loadings, Y_scores, data_loadings, as.matrix(data_scores)), "data_scores is not a tibble")
})