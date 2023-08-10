
test_that("input: mat must be matrix", {
  expect_error(normalize_loadings(mat = as.data.frame(matrix(1:9, 3, 3)), type = "largest_abs"), "mat is not a matrix.")
})

test_that("input: type invalid", {
  expect_error(normalize_loadings(mat = matrix(1:9, 3, 3), type = "something else"), "type specified is not valid.")
})

test_that("functionality: with type largest_abs", {
  expect_equal(normalize_loadings(mat = cbind(c(-1, 0, 1), c(-2, 0, 2), c(-3, 0, 3)), type = "largest_abs"), matrix(rep(c(-1, 0, 1), 3), 3, 3))
  expect_equal(normalize_loadings(mat = cbind(c(-2, 0, 1), c(2, 0, 1)), type = "largest_abs"), cbind(c(-1, 0, 0.5), c(1, 0, 0.5)) )
})

