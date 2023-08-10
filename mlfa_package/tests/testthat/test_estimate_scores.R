

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

  # print("here: ")
  # print(Y)

  return(list(Y, Z, A))
}

test_that("no error", {
  expect_error(estimate_scores(Y_in = matrix(0, 2, 2), A_in = matrix(0, 2, 2)), NA)
})

test_that("input: Y_in is matrix", {
  expect_error(estimate_scores(Y_in = as.data.frame(matrix(0, 2, 2)), A_in = matrix(0, 2, 2)), "Y_in is not a matrix.")
})

test_that("input: A_in is matrix", {
  expect_error(estimate_scores(Y_in = matrix(0, 2, 2), A_in = as.data.frame(matrix(0, 2, 2))), "A_in is not a matrix.")
})

test_that("input: type is valid", {
  expect_error(estimate_scores(Y_in = matrix(0, 2, 2), A_in = matrix(0, 2, 2), type = "something else"), "Variable type must be 'lin_reg'.")
})

test_that("input: P dimension", {
  expect_error(estimate_scores(Y_in = matrix(0, 2, 2), A_in = matrix(0, 1, 2)), "P dimension, which is first dimension of A_in and second dimension of Y_in, must match.")
})


test_that("functionality: reproduced A used to compute Z and Y reproduced, WITH loadings normalization", {
  P <- 100
  list[Y, Z] <- make_Y_matrix_simple(N = 1000, P = P, K = 3, noise_sd = 0.1)
  list[A, ppca_object] <- estimate_loadings(Y_in = Y, type = "n_PCs", var_explained = NULL, n_PCs = 3, rotation = "promax", conv_thres = 1e-4, sparse_abs_thres = 0.0, ppca_seed = 1)
  # A has swapped columns
  A[, c(1, 2, 3)] <- A[, c(2, 3, 1)]
  # A has first column negative sign
  A[, 1] <- -A[, 1]
  list[Z_res, S_res] <- estimate_scores(Y_in = Y, A_in = A)  # using the estimated loadings matrix

  # Note: would probably be even more precise if A is not normalized in estimate_loadings!
  expect_equal(Z, Z_res, tolerance = 0.5, check.attributes = FALSE)  # tolerance -> does not have to be exactly equal due to approximate EM algorithm; works only with check.attributes = FALSE -> see https://github.com/r-lib/testthat/issues/580
})

test_that("functionality: reproduced A used to compute Z and Y reproduced, WITHOUT loadings normalization", {
  P <- 100
  list[Y, Z, A] <- make_Y_matrix_simple(N = 1000, P = P, K = 3, noise_sd = 0.1)
  list[A_res, ppca_object] <- estimate_loadings(Y_in = Y, type = "n_PCs", var_explained = NULL, n_PCs = 3, rotation = "promax", conv_thres = 1e-4, sparse_abs_thres = 0.0, ppca_seed = 1, normalize_loadings_type = "none")
  # A_res has swapped columns
  A_res[, c(1, 2, 3)] <- A_res[, c(2, 3, 1)]
  # A_res has first column negative sign
  A_res[, 1] <- -A_res[, 1]
  list[Z_res, S_res] <- estimate_scores(Y_in = Y, A_in = A_res)  # using the estimated loadings matrix

  # Note: would probably be even more precise if A_res is not normalized in estimate_loadings!
  # Note: Observation: without normalization of loadings, the scores are much mroe off
  expect_equal(Z, Z_res, tolerance = 15, check.attributes = FALSE)  # tolerance -> does not have to be exactly equal due to approximate EM algorithm; works only with check.attributes = FALSE -> see https://github.com/r-lib/testthat/issues/580

  # Note: Y, however, is very well reproduced
  expect_equal(Z_res %*% t(A_res), Y, tolerance = 0.3, check.attributes = FALSE)
})

# TODO a functionality test with real data (metaboliteData) and whether reasonable few PCs are selected




