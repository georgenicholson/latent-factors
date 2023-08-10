

test_that("input: data is a tibble", {
  expect_error(as_ranks(data = data.frame(a = 1:10, b = 11:20, c=31:40), columns = c("b", "c")))
})

test_that("input: columns are in data", {
  expect_error(as_ranks(data = dplyr::tibble(a = 1:10, b = 11:20, c=31:40), columns = c("d")))
})

test_that("functionality: decimal ranks correctly computed", {
  expect_equal(as_ranks(data = dplyr::tibble(a = 1:10, b = 11:20, c=31:40), columns = c("b", "c")),
               dplyr::tibble(a = 1:10, b = round((1:10)/10-0.1, 1), c = round((1:10)/10-0.1, 1) ))  # rounding due to precision issues
})

test_that("functionality: works with negative values equally", {
  expect_equal(as_ranks(data = dplyr::tibble(a = 1:10, b = -9:0, c=31:40), columns = c("b", "c")),
               dplyr::tibble(a = 1:10, b = round((1:10)/10-0.1, 1), c = round((1:10)/10-0.1, 1) ))  # rounding due to precision issues
})

test_that("functionality: ties are averaged", {
  expect_equal(as_ranks(data = dplyr::tibble(a = 1:4, b = c(1, 2, 2, 3)), columns = "b"),
               dplyr::tibble(a = 1:4, b = c(0, 0.375, 0.375, 0.75)) )
})

test_that("functionality: NAs are ignored", {
  expect_identical(as_ranks(data = dplyr::tibble(a = 1:5, b = c(1, 2, NA, 3, 4)), columns = "b"),
               dplyr::tibble(a = 1:5, b = c(0.0, 0.25, NA, 0.5, 0.75)))  # tolerance allows to ignore small precision issues -> use also above
})



