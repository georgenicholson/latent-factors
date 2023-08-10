

test_that("input: mat is matrix", {
  expect_error(sparse_heur_thres(mat = as.data.frame(matrix(c(0.1, 0.2, 0.3, 0.4), 2, 2)), thres = 0.1), "mat is not a matrix.")
})


test_that("input: thres must be non-negative.", {
  expect_error(sparse_heur_thres(mat = matrix(c(0.1, 0.2, 0.3, 0.4), 2, 2), thres = -0.1), "thres must be non-negative double.")
})

test_that("functionality: correctly sets values to 0 (also checking edge case)", {
  expect_equal(sparse_heur_thres(mat = matrix(c(0.1, 0.2, 0.3, 0.4), 2, 2), thres = 0.1), matrix(c(0.1, 0.2, 0.3, 0.4), 2, 2))
  expect_equal(sparse_heur_thres(mat = matrix(c(0.1, 0.2, 0.3, 0.4), 2, 2), thres = 0.11), matrix(c(0, 0.2, 0.3, 0.4), 2, 2))
  expect_equal(sparse_heur_thres(mat = matrix(c(-0.1, 0.2, 0.3, 0.4), 2, 2), thres = 0.11), matrix(c(0, 0.2, 0.3, 0.4), 2, 2))
  expect_equal(sparse_heur_thres(mat = matrix(c(-0.1, 0.2, 0.3, 0.4), 2, 2), thres = 0.1), matrix(c(-0.1, 0.2, 0.3, 0.4), 2, 2))
})


