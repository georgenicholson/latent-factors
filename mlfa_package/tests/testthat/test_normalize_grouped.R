test_that("input: data must be tibble", {
  expect_error(normalize_grouped(data = data.frame(a = c(1, NA, 3), b = c(1, NA, 3), c = c(1, 2, 3)), columns = c("a", "c"), type = "standardize", group = NULL), "data is not a tibble.")
})


test_that("input: columns valid", {
  expect_error(normalize_grouped(data = tibble(a = c(1, NA, 3), b = c(1, NA, 3), c = c(1, 2, 3)), columns = c("a", "d"), type = "standardize", group = NULL))
})

test_that("input: type valid", {
  expect_error(normalize_grouped(data = tibble(a = c(1, NA, 3), b = c(1, NA, 3), c = c(1, 2, 3)), columns = c("a", "c"), type = "something else", group = NULL))
})

test_that("input: group valid", {
  expect_error(normalize_grouped(data = tibble(a = c(1, NA, 3), b = c(1, NA, 3), c = c(1, 2, 3)), columns = c("a", "c"), type = "standardize", group = "d"))
})

test_that("functionality: correct normalization (without group, standardize and center)", {
  # standardize
  data_norm <- normalize_grouped(data = tibble(a = c(1, NA, 3), b = c(1, NA, 3), c = c(0, 2, 4)), columns = c("c"), type = "standardize", group = NULL)
  expect_equal(data_norm, tibble(a = c(1, NA, 3), b = c(1, NA, 3), c = c(-1, 0, 1)))

  data_norm <- normalize_grouped(data = tibble(a = c(1, 2, 3), b = c(1, NA, 3), c = c(0, 2, 4)), columns = c("a"), type = "standardize", group = NULL)
  expect_equal(data_norm, tibble(a = c((1-mean(c(1, 2, 3)))/sd(c(1, 2, 3)), (2-mean(c(1, 2, 3)))/sd(c(1, 2, 3)), (3-mean(c(1, 2, 3)))/sd(c(1, 2, 3))), b = c(1, NA, 3), c = c(0, 2, 4)))  # x - mu / sd , sd is unbiased (N-1 in denominator)

  # center
  data_norm <- normalize_grouped(data = tibble(a = c(1, NA, 3), b = c(1, NA, 3), c = c(0, 2, 4)), columns = c("c"), type = "center", group = NULL)
  expect_equal(data_norm, tibble(a = c(1, NA, 3), b = c(1, NA, 3), c = c(-2, 0, 2)))
})


test_that("functionality: correct normalization (with group, standardize)", {
  # center
  data_norm <- normalize_grouped(data = tibble(a = c(1, NA, 3, NA, NA, NA), b = c(1, 2, 3, 1, 2, 3), c = c("gr1", "gr1", "gr1", "gr2", "gr2", "gr2")), columns = c("b"), type = "standardize", group = c("c"))
  expect_equal(data_norm, tibble(a = c(1, NA, 3, NA, NA, NA), b = c(-1, 0, 1, -1, 0, 1), c = c("gr1", "gr1", "gr1", "gr2", "gr2", "gr2")) )
})




# TODO test if multiple columns possible
# TODO test if multiple variables in group possible
