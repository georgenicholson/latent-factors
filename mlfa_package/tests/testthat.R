# must first load package with: 
# devtools::load_all("../../mlfa_package")

# to test the package: 
# set working directory to package
# then: devtools::test()

library(testthat)
library(mlfa)

test_check("mlfa")
