# loads the data into a variable called "metaboliteData"
utils::data(metaboliteData)
# transform the "metaboliteData" object into a tibble
data <- dplyr::tibble(as.data.frame(metaboliteData))

# var_group_1 <- colnames(metaboliteData) # [1:20]
# var_group_2 <- colnames(metaboliteData) # [21:52]
# var_groups <- list(var_group_1, var_group_2)

# must add
data$subj_ID <- rep('subj_1', dim(data)[1])
data$visit_ID <- rep(1, dim(data)[1])

vars = colnames(metaboliteData)

# extra char column
data$char_col <- rep("a", dim(data)[1])

test_that("no errors", {
  expect_error(mfd_object <- mfd(data = data, vars = vars, time = 'visit_ID', subject = 'subj_ID', loadings.n_PCs = 3),
               NA)
})

test_that("input: data must be tibble", {
  expect_error(mfd_object <- mfd(data = as.data.frame(data), vars = vars, time = 'visit_ID', subject = 'subj_ID', loadings.n_PCs = 3),
               "data is not a tibble.")
})

test_that("input: Less than two columns present in data -> cannot reasonable estimate PCs.", {
  expect_error(mfd_object <- mfd(data = data[, 1], vars = vars, time = 'visit_ID', subject = 'subj_ID'),
               "Less than two columns present in data -> cannot reasonable estimate PCs.")
})

test_that("input: vars must be columns of data", {
  vars_diff = vars
  vars_diff[1] = "something else"
  expect_error(mfd_object <- mfd(data = data, vars = vars_diff, time = 'visit_ID', subject = 'subj_ID'),
               "vars must be columns of data.")
})


test_that("input: vars columns must be numeric", {
  expect_error(mfd_object <- mfd(data = data, vars = c(vars, "char_col"), time = 'visit_ID', subject = 'subj_ID', loadings.n_PCs = 3),
               "char_col is not numeric.")
})

test_that("input: vars must be columns of data", {
  vars_diff = vars
  vars_diff[1] = "something else"
  expect_error(mfd_object <- mfd(data = data, vars = vars_diff, time = 'visit_ID', subject = 'subj_ID'),
               "vars must be columns of data.")
})

test_that("input: time must be columns of data", {
  expect_error(mfd_object <- mfd(data = data, vars = vars, time = 'visit_ID', subject = 'subj_ID'),
               NA)
  expect_error(mfd_object <- mfd(data = data, vars = vars, time = 'something else', subject = 'subj_ID'),
               "time must be column of data.")
})

test_that("input: subject must be columns of data", {
  expect_error(mfd_object <- mfd(data = data, vars = vars, time = 'visit_ID', subject = 'subj_ID'),
               NA)
  expect_error(mfd_object <- mfd(data = data, vars = vars, time = 'visit_ID', subject = 'something else'),
               "subject must be column of data.")
})


test_that("inputs: loadings.n_PCs with n_PCs", {
  expect_error(mfd_object <- mfd(data = data, vars = vars, time = 'visit_ID', subject = 'subj_ID', loadings.n_PCs = "max"),
               NA)
  expect_error(mfd_object <- mfd(data = data, vars = vars, time = 'visit_ID', subject = 'subj_ID', loadings.n_PCs = 1),
               "loadings.n_PCs must be either 'max' or a positive integer >= 2 and < length of vars.")
  expect_error(mfd_object <- mfd(data = data, vars = vars, time = 'visit_ID', subject = 'subj_ID', loadings.n_PCs = 2.5),
               "loadings.n_PCs must be either 'max' or a positive integer >= 2 and < length of vars.")
  expect_error(mfd_object <- mfd(data = data, vars = vars, time = 'visit_ID', subject = 'subj_ID', loadings.n_PCs = -1),
               "loadings.n_PCs must be either 'max' or a positive integer >= 2 and < length of vars.")
  expect_error(mfd_object <- mfd(data = data, vars = vars, time = 'visit_ID', subject = 'subj_ID', loadings.n_PCs = length(vars)),
               "loadings.n_PCs must be either 'max' or a positive integer >= 2 and < length of vars.")
  expect_error(mfd_object <- mfd(data = data, vars = vars, time = 'visit_ID', subject = 'subj_ID', loadings.n_PCs = 2),
               NA)
})


test_that("input: loadings.sparse_abs_thres within 0.0 and 1.0", {
  # edge cases
  expect_error(mfd_object <- mfd(data = data, vars = vars, time = 'visit_ID', subject = 'subj_ID', loadings.n_PCs = 3, loadings.sparse_abs_thres = 0.0),
               NA)
  expect_error(mfd_object <- mfd(data = data, vars = vars, time = 'visit_ID', subject = 'subj_ID', loadings.n_PCs = 3, loadings.sparse_abs_thres = 0.9),
               NA)
  expect_error(mfd_object <- mfd(data = data, vars = vars, time = 'visit_ID', subject = 'subj_ID', loadings.n_PCs = 3, loadings.sparse_abs_thres = 1.5),
               "sparse_abs_thres must be double between 0.0 and 1.0, including bounds.")
  expect_error(mfd_object <- mfd(data = data, vars = vars, time = 'visit_ID', subject = 'subj_ID', loadings.n_PCs = 3, loadings.sparse_abs_thres = -1.5),
               "sparse_abs_thres must be double between 0.0 and 1.0, including bounds.")
})


test_that("input: seed correctly specified", {
  expect_error(mfd_object <- mfd(data = data, vars = vars, time = 'visit_ID', subject = 'subj_ID', loadings.n_PCs = 3, seed = NULL),
               NA)
  expect_error(mfd_object <- mfd(data = data, vars = vars, time = 'visit_ID', subject = 'subj_ID', loadings.n_PCs = 3, seed = 1.5),
               "seed must be NULL or positive interger.")
  expect_error(mfd_object <- mfd(data = data, vars = vars, time = 'visit_ID', subject = 'subj_ID', loadings.n_PCs = 3, seed = -1),
               "seed must be NULL or positive interger.")
})

test_that("input: loadings.rotation correctly specified", {
  expect_error(mfd_object <- mfd(data = data, vars = vars, time = 'visit_ID', subject = 'subj_ID', loadings.n_PCs = 3, loadings.rotation = "promax"),
               NA)
  expect_error(mfd_object <- mfd(data = data, vars = vars, time = 'visit_ID', subject = 'subj_ID', loadings.n_PCs = 3, loadings.rotation = "varimax"),
               NA)
  expect_error(mfd_object <- mfd(data = data, vars = vars, time = 'visit_ID', subject = 'subj_ID', loadings.rotation = "something else"),
               "loadings.rotation must be either 'varimax' or 'promax'.")
})

# TODO functionality and all attribute check!

