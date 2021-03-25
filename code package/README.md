# what is this?
the folder structure that of an R-package. to compile the package in Rstudio run *[Build]>[Clean and rebuild]*. afterwards all (read: both) functions are avialable via `library(riiv)`

# which files are important?

1. scripts/example simulation.R
containts the code to conduct the example simulations that uses  2 and 3

2. src/ri_iv_test_cpp.cpp
defines the function to compute randomization inference p-values

3. R/ri_iv_gridsearch.R
defines the function for the gridsearch to estimate RI-based confidence bands/sets
