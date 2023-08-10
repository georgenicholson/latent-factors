
test_A <- matrix(1:9, 3, 3)
test_Z <- matrix(1:9, 3, 3)
test_S <- matrix(1:9, 3, 3)

test_that("no error", {
  expect_error(normalize_stage_1_2_results(A_in = test_A, Z_in = test_Z, S_in = test_S, type = "normalize_custom"), NA)
  expect_error(normalize_stage_1_2_results(A_in = test_A, Z_in = test_Z, S_in = test_S, type = "scale"), NA)
})

test_that("input: A_in is matrix", {
  expect_error(normalize_stage_1_2_results(A_in = as.data.frame(test_A), Z_in = test_Z, S_in = test_S, type = "normalize_custom"), "A_in is not a matrix.")
})

test_that("input: Z_in is matrix", {
  expect_error(normalize_stage_1_2_results(A_in = test_A, Z_in = as.data.frame(test_Z), S_in = test_S, type = "normalize_custom"), "Z_in is not a matrix.")
})

test_that("input: S_in is matrix", {
  expect_error(normalize_stage_1_2_results(A_in = test_A, Z_in = test_Z, S_in = as.data.frame(test_S), type = "normalize_custom"), "S_in is not a matrix.")
})


test_that("no error", {
  expect_error(normalize_stage_1_2_results(A_in = test_A, Z_in = test_Z, S_in = test_S, type = "something else"), "type specifies no valid normalization.")
})


# TODO functionality test
