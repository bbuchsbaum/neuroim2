library(testthat)

# Test cases for split_reduce
test_that("split_reduce works with matrix, integer and function inputs", {
  # Create input matrix
  mat <- matrix(1:12, nrow=4)

  # Split matrix by factor
  fac <- as.factor(c(1, 1, 2, 2))

  # Define custom function
  custom_fun <- function(x) sum(x)

  # Test split_reduce
  output <- split_reduce(mat, fac, custom_fun)

  # Expected output
  expected_output <- matrix(c(3, 7, 11, 15,19,23), nrow=2)

  expect_equal(output, expected_output)
})

# Test cases for split_scale
test_that("split_scale works with matrix and factor inputs", {
  # Create input matrix
  mat <- matrix(1:12, nrow=4)

  # Split matrix by factor
  fac <- as.factor(c(1, 1, 2, 2))

  # Test split_scale
  output <- split_scale(mat, fac, center=TRUE, scale=TRUE)

  # Expected output
  expected_output <- rbind(scale(mat[1:2,], center=TRUE, scale=TRUE),
                           scale(mat[3:4,], center=TRUE, scale=TRUE))

  expect_equal(output, expected_output)
})




