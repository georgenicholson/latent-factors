
# input_test_Y <- matrix(0, 100, 10)
# input_test_A <- matrix(0, 3, 10)
#
# # NA -> this case means we expext no error for this call
# # TODO implement such a test case for every single test file
# test_that("no error", {
#   expect_error(estimate_Z_S_helper(Y_in = input_test_Y, A_in = input_test_A), NA)
# })
#
# test_that("input: Y_in is matrix", {
#           expect_error(estimate_Z_S_helper(Y_in = data.frame(input_test_Y), A_in = input_test_A), "Y_in is not a matrix.")
#           })
#
# test_that("input: A_in is matrix", {
#   expect_error(estimate_Z_S_helper(Y_in = input_test_Y, A_in = data.frame(input_test_A)), "A_in is not a matrix.")
# })
#
# test_that("input: A_in is matrix", {
#   expect_error(estimate_Z_S_helper(Y_in = input_test_Y, A_in = matrix(0, 3, 9)), "Second dimension of Y_in and A_in must be same.")
# })


# TODO functionality tests
