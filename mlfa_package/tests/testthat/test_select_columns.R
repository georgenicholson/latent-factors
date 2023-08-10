test_that("input: data is tibble", {
  expect_error(select_columns(data = data.frame(a = 1:3, b = 1:3, c = 1:3), columns = c("c")) )
})

test_that("input: selected_columns in data", {
  expect_error(select_columns(data = tibble(a = 1:3, b = 1:3, c = 1:3), columns = c("d")) )
})

test_that("functionality: correct columns selected, can select multiple columns", {
  expect_equal(select_columns(data = tibble(a = 1:3, b = 1:3, c = 1:3), columns = c("a", "c")), tibble(a = 1:3, c = 1:3))
})
