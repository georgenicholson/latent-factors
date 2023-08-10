
Y_test = rbind(c(NA, NA, 1),
               c(NA, NA, 2.5),
               c(NA, 1.2, 2.4))

test_that("no error", {
  expect_error(find_unique_missing_patterns(Y_in = Y_test), NA)
})

test_that("input: Y_in is matrix", {
  expect_error(find_unique_missing_patterns(Y_in = data.frame(Y_test)), "Y_in is not a matrix.")
})


test_that("functionality: both outputs as expected", {
  list[unique_missing_patterns_char, indices_Y_in_list] <- find_unique_missing_patterns(Y_test)
  expect_equal(unique_missing_patterns_char, c("100", "110"))
  expect_equal(indices_Y_in_list, list(c("3"), c("1", "2")))
})


