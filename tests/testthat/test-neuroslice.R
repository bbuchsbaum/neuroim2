library(testthat)
library(neuroim2)

# Test the NeuroSlice constructor
test_that("NeuroSlice constructor works correctly", {
  data <- matrix(rnorm(64 * 64), 64, 64)
  space <- NeuroSpace(c(64, 64), spacing = c(1, 1))
  nslice <- NeuroSlice(data, space)

  expect_s4_class(nslice, "NeuroSlice")
  expect_equal(nslice@space, space)
  expect_equal(nslice@.Data, data)
})

# Test grid_to_index methods
test_that("grid_to_index methods work correctly", {
  data <- matrix(rnorm(64 * 64), 64, 64)
  space <- NeuroSpace(c(64, 64), spacing = c(1, 1))
  nslice <- NeuroSlice(data, space)

  coords <- matrix(c(1, 1, 64, 64), nrow=2, byrow=TRUE)
  idx <- grid_to_index(nslice, coords)

  expect_equal(idx, c(1, 4096))

  coords <- c(1, 1)
  idx <- grid_to_index(nslice, coords)

  expect_equal(idx, 1)
})

# Test index_to_grid method
test_that("index_to_grid method works correctly", {
  data <- matrix(rnorm(64 * 64), 64, 64)
  space <- NeuroSpace(c(64, 64), spacing = c(1, 1))
  nslice <- NeuroSlice(data, space)

  idx <- c(1, 4096)
  coords <- index_to_grid(nslice, idx)

  expect_equal(coords, matrix(c(1, 1, 64, 64), nrow = 2, byrow = TRUE))
})


# Test the mapToColors function
test_that("mapToColors works correctly", {
  data <- matrix(rnorm(64 * 64), 64, 64)
  data[1] = 0
  col <- heat.colors(128, alpha = 1)
  zero_col <- "#00000000"
  alpha <- 1
  irange <- range(data)
  threshold <- c(-2, 2)

  mapped_colors <- mapToColors(data, col, zero_col, alpha=1, irange)

  expect_equal(dim(mapped_colors), dim(data))
  expect_equal(mapped_colors[data == 0], zero_col)

  threshold <- c(-1, 1)
  mapped_colors <- mapToColors(data, col, zero_col, alpha=1, irange, threshold)

  expect_equal(dim(mapped_colors), dim(data))
  expect_equal(mapped_colors[data == 0], zero_col)

  n <- sum(data >= threshold[1] & data <= threshold[2])
  expect_equivalent(mapped_colors[data >= threshold[1] & data <= threshold[2]], rep("#00000000", n))
})


