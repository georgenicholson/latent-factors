
test_that("input: mat is matrix", {
  expect_error(sparse_rotation(mat = as.data.frame(matrix(1:4, 2, 2)), type = "promax"), "mat is not a matrix.")
})

test_that("input: mat is matrix", {
  expect_error(sparse_rotation(mat = matrix(1:4, 2, 2), type = "something else"), "type is not a valid rotation type.")
})

test_that("input: mat is matrix", {
  expect_error(sparse_rotation(mat = matrix(1:4, 4, 1), type = "promax"), "mat has less than two columns so that rotation cannot be performed.")
})


# TODO functionality test
