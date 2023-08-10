make_Y_matrix_complex <- function(N, P, K, noise_sd) {
  set.seed(5)
  A <- rbind(c(rep(1, 10), rep(0, 90)),
             c(rep(0, 10), rep(1, 10), rep(0, 80)),
             c(rep(0, 20), rep(1, 80))
             )
  # print(dim(A))
  # assert(dim(A) == c(K, P))
  Z <- NULL
  for (i in 1:N) {
    new_row <- seq(from = 0, to = 0.01, length = K) * i
    Z <- rbind(Z, new_row)
  }
  # Z <- matrix(rnorm(N*K), nrow = N, ncol = K)
  # print(dim(Z))
  # assert(dim(Z) == c(N, K))
  E <- rnorm(n = c(N, P), mean = 0, sd = noise_sd)
  Y <- Z %*% A + E

  return(Y)
}

#'   Y   =   Z   %*%  t(A)  +   E
#' (NxP) = (NxK) %*% (KxP)  + (NxP)
#' Note: E is implicit in add Gaussian noise call.
make_Y_matrix_simple <- function(N = 1000, P = 100, K = 3, noise_sd = 0.1) {
  set.seed(5)
  # A is PxK
  # non-zero values ensures A is not singular
  A <- cbind(
    c(1, rep(0.1, P-1)),
    c(0.05, 1, rep(0.2, P-2)),
    c(0.1, 0.1, rep(1, P-2))
  )
  # print(paste0("dim(A): ", dim(A)))
  # assert(dim(A) == c(K, P))
  Z <- matrix(rnorm(n = N*K, mean = 0), nrow = N, ncol = K)
  # print(paste0("dim(Z): ", dim(Z)))
  # assert(dim(Z) == c(N, K))
  E <- matrix(rnorm(n = N*P, mean = 0, sd = noise_sd), nrow = N, ncol = P)
  Y <- Z %*% t(A) + E

  return(Y)
}

input_testing_Y <- matrix(rnorm(1000*100), 1000, 100)

test_that("no error", {
  # with toy matrix
  expect_error(
    estimate_loadings(Y_in = input_testing_Y, type = "var_explained", var_explained = 0.5, n_PCs = NULL, rotation = "promax", conv_thres = 1e-4, sparse_abs_thres = +1.0, ppca_seed = 1),
    NA
  )
  # with real data
  # loads the data into a variable called "metaboliteData"
  data(metaboliteData)
  # transform the "metaboliteData" object into a matrix
  data <- as.matrix(dplyr::tibble(as.data.frame(metaboliteData)))
  # throws warning (regarding precision), but no error. Warning message: "Precision for components 50 - 51 is below .Machine$double.eps. Results for those components are likely to be inaccurate!!"
  # we surpress this warning for now. -> before: expect_warning(...); right now, this unit test is not doing anything.
  list[A, ppca_object] <- estimate_loadings(Y_in = data, type = "n_PCs", var_explained = NULL, n_PCs = 3, rotation = "promax", conv_thres = 1e-4, sparse_abs_thres = 0.0, ppca_seed = 1)
})

test_that("input: Y_in is matrix", {
  expect_error(
    estimate_loadings(Y_in = as.data.frame(input_testing_Y), type = "var_explained", var_explained = 0.5, n_PCs = NULL, rotation = "promax", conv_thres = 1e-4, sparse_abs_thres = 0.0, ppca_seed = 1),
    "Y_in is not a matrix."
  )
})

test_that("input: type", {
  expect_error(
    estimate_loadings(Y_in = input_testing_Y, type = "something else", var_explained = 0.5, n_PCs = NULL, rotation = "promax", conv_thres = 1e-4, sparse_abs_thres = 0.0, ppca_seed = 1),
    "type must be either 'var_explained' or 'n_PCs'."
  )
})

test_that("input: var_explained", {
  expect_error(
    estimate_loadings(Y_in = input_testing_Y, type = "var_explained", var_explained = 1.5, n_PCs = NULL, rotation = "promax", conv_thres = 1e-4, sparse_abs_thres = 0.0, ppca_seed = 1),
    "var_explained must be double in > 0.0 and <= 1.0."
  )
  expect_error(
    estimate_loadings(Y_in = input_testing_Y, type = "var_explained", var_explained = -0.5, n_PCs = NULL, rotation = "promax", conv_thres = 1e-4, sparse_abs_thres = 0.0, ppca_seed = 1),
    "var_explained must be double in > 0.0 and <= 1.0."
  )
})

test_that("input: n_PCs", {
  expect_error(
    estimate_loadings(Y_in = input_testing_Y, type = "n_PCs", var_explained = NULL, n_PCs = 2.5, rotation = "promax", conv_thres = 1e-4, sparse_abs_thres = 0.0, ppca_seed = 1),
    "n_PCs must be positive integer >= 2 and < P."
  )
  expect_error(
    estimate_loadings(Y_in = input_testing_Y, type = "n_PCs", var_explained = NULL, n_PCs = -2, rotation = "promax", conv_thres = 1e-4, sparse_abs_thres = 0.0, ppca_seed = 1),
    "n_PCs must be positive integer >= 2 and < P."
  )
  expect_error(
    estimate_loadings(Y_in = input_testing_Y, type = "n_PCs", var_explained = NULL, n_PCs = dim(input_testing_Y)[2], rotation = "promax", conv_thres = 1e-4, sparse_abs_thres = 0.0, ppca_seed = 1),
    "n_PCs must be positive integer >= 2 and < P."
  )
})

test_that("input: type and parameter specifying combination", {
  expect_error(
    estimate_loadings(Y_in = input_testing_Y, type = "var_explained", var_explained = 0.5, n_PCs = 2, rotation = "promax", conv_thres = 1e-4, sparse_abs_thres = 0.0, ppca_seed = 1),
    "type is 'var_explained', i.e. only var_explained should be specified, n_PCs should be NULL."
  )
  expect_error(
    estimate_loadings(Y_in = input_testing_Y, type = "var_explained", var_explained = NULL, n_PCs = 2, rotation = "promax", conv_thres = 1e-4, sparse_abs_thres = 0.0, ppca_seed = 1),
    "type is 'var_explained', i.e. only var_explained should be specified, n_PCs should be NULL."
  )
  expect_error(
    estimate_loadings(Y_in = input_testing_Y, type = "n_PCs", var_explained = 0.5, n_PCs = 2, rotation = "promax", conv_thres = 1e-4, sparse_abs_thres = 0.0, ppca_seed = 1),
    "type is 'n_PCs', i.e. only n_PCs should be specified, var_explained should be NULL."
  )
  expect_error(
    estimate_loadings(Y_in = input_testing_Y, type = "n_PCs", var_explained = 0.5, n_PCs = NULL, rotation = "promax", conv_thres = 1e-4, sparse_abs_thres = 0.0, ppca_seed = 1),
    "type is 'n_PCs', i.e. only n_PCs should be specified, var_explained should be NULL."
  )
})

test_that("input: rotation", {
  expect_error(
    estimate_loadings(Y_in = input_testing_Y, type = "var_explained", var_explained = 0.5, n_PCs = NULL, rotation = "something else", conv_thres = 1e-4, sparse_abs_thres = 0.0, ppca_seed = 1),
    "rotation must be either 'varimax' or 'promax'."
  )
})

test_that("input: conv_thres", {
  expect_error(
    estimate_loadings(Y_in = input_testing_Y, type = "var_explained", var_explained = 0.5, n_PCs = NULL, rotation = "promax", conv_thres = -1e-4, sparse_abs_thres = 0.0, ppca_seed = 1),
    "conv_thres must be positive double."
  )
})

test_that("input: sparse_thres", {
  expect_error(
    estimate_loadings(Y_in = input_testing_Y, type = "var_explained", var_explained = 0.5, n_PCs = NULL, rotation = "promax", conv_thres = 1e-4, sparse_abs_thres = -0.1, ppca_seed = 1),
    "sparse_abs_thres must be double between 0.0 and 1.0, including bounds."
  )
  expect_error(
    estimate_loadings(Y_in = input_testing_Y, type = "var_explained", var_explained = 0.5, n_PCs = NULL, rotation = "promax", conv_thres = 1e-4, sparse_abs_thres = 1.1, ppca_seed = 1),
    "sparse_abs_thres must be double between 0.0 and 1.0, including bounds."
  )
  expect_error(
    estimate_loadings(Y_in = input_testing_Y, type = "var_explained", var_explained = 0.5, n_PCs = NULL, rotation = "promax", conv_thres = 1e-4, sparse_abs_thres = 2.1, ppca_seed = 1),
    "sparse_abs_thres must be double between 0.0 and 1.0, including bounds."
  )
  # 0.0 is ok
  expect_error(
    estimate_loadings(Y_in = input_testing_Y, type = "var_explained", var_explained = 0.5, n_PCs = NULL, rotation = "promax", conv_thres = 1e-4, sparse_abs_thres = 0.0, ppca_seed = 1),
    NA
  )
})

test_that("input: ppca_seed", {
  expect_error(
    estimate_loadings(Y_in = input_testing_Y, type = "var_explained", var_explained = 0.5, n_PCs = NULL, rotation = "promax", conv_thres = 1e-4, sparse_abs_thres = 0.1, ppca_seed = 1.5),
    "ppca_seed must be positive integer >= 1."
  )
  expect_error(
    estimate_loadings(Y_in = input_testing_Y, type = "var_explained", var_explained = 0.5, n_PCs = NULL, rotation = "promax", conv_thres = 1e-4, sparse_abs_thres = 0.1, ppca_seed = -5),
    "ppca_seed must be positive integer >= 1."
  )
})

test_that("input: Y_in at least 2 variables", {
  expect_error(
    estimate_loadings(Y_in = as.matrix(input_testing_Y[, 1]), type = "var_explained", var_explained = 0.5, n_PCs = NULL, rotation = "promax", conv_thres = 1e-4, sparse_abs_thres = 0.1, ppca_seed = 1),
    "Less than two columns present in Y_in -> cannot reasonably estimate PCs."
  )
})

test_that("input: n_PCs < P (already covered above)", {
  expect_error(
    estimate_loadings(Y_in = input_testing_Y, type = "n_PCs", var_explained = NULL, n_PCs = "max", rotation = "promax", conv_thres = 1e-4, sparse_abs_thres = 0.1, ppca_seed = 1),
    NA
  )
  expect_error(
    estimate_loadings(Y_in = input_testing_Y[, c(1,2)], type = "n_PCs", var_explained = NULL, n_PCs = 3, rotation = "promax", conv_thres = 1e-4, sparse_abs_thres = 0.1, ppca_seed = 1),
    "n_PCs must be positive integer >= 2 and < P or 'max'."
  )
  expect_error(
    estimate_loadings(Y_in = input_testing_Y[, c(1,2)], type = "n_PCs", var_explained = NULL, n_PCs = 1, rotation = "promax", conv_thres = 1e-4, sparse_abs_thres = 0.1, ppca_seed = 1),
    "n_PCs must be positive integer >= 2 and < P or 'max'."
  )
})

test_that("input: valid normalize_loadings_type", {
  expect_error(
    estimate_loadings(Y_in = input_testing_Y, type = "n_PCs", var_explained = NULL, n_PCs = 3, rotation = "promax", conv_thres = 1e-4, sparse_abs_thres = 0.1, ppca_seed = 1, normalize_loadings_type = "something else"),
    "normalize_loadings_type is not valid."
  )
})





# IMPORTANT NOTE: estimate_loadings might throw errors like
# "Error in if (rel_ch < threshold & count > 5) { :
#    missing value where TRUE/FALSE needed"
# This is because nPCs in the ppca call is set to high
# TODO how to resolve this? -> alwyas decrease nPCs consecutively?
test_that("functionality: simulated A with expected dimensions can be reproduced", {
  P <- 100
  list[A, ppca_object] <- estimate_loadings(Y_in = make_Y_matrix_simple(N = 1000, P = P, K = 3, noise_sd = 0.1), type = "n_PCs", var_explained = NULL, n_PCs = 3, rotation = "promax", conv_thres = 1e-4, sparse_abs_thres = 0.0, ppca_seed = 1)
  # A has swapped columns
  A[, c(1, 2, 3)] <- A[, c(2, 3, 1)]
  # A has first column negative sign
  A[, 1] <- -A[, 1]
  A_exp <- cbind(
    c(1, rep(0.1, P-1)),
    c(0.05, 1, rep(0.2, P-2)),
    c(0.1, 0.1, rep(1, P-2)))
  expect_equal(A, A_exp, tolerance = 0.3, check.attributes = FALSE)  # tolerance -> does not have to be exactly equal due to approximate EM algorithm; works only with check.attributes = FALSE -> see https://github.com/r-lib/testthat/issues/580
})

# TODO functionality tests

# TODO unit test functionality of rotation
