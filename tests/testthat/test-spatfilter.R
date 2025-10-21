library(testthat)
library(neuroim2)

# Test the gaussian_blur function
test_that("gaussian_blur works correctly", {
  vol <- NeuroVol(array(rnorm(64 * 64 * 64), c(64, 64, 64)), NeuroSpace(c(64, 64, 64)))
  mask <- as.logical(NeuroVol(as.logical(round(runif(prod(dim(vol))))), space(vol)))
  sigma <- 2
  window <- 1

  blurred_vol <- gaussian_blur(vol, mask, sigma, window)

  expect_s4_class(blurred_vol, "NeuroVol")
  expect_equal(dim(blurred_vol), dim(vol))
  expect_true(all(blurred_vol[mask == 0] == 0))
})

# Test the guided_filter function
test_that("guided_filter works correctly", {
  vol <- NeuroVol(array(rnorm(32 * 32 * 32), c(32, 32, 32)), NeuroSpace(c(32, 32, 32)))
  radius <- 4
  epsilon <- 0.7^2

  filtered_vol <- guided_filter(vol, radius, epsilon)

  expect_s4_class(filtered_vol, "NeuroVol")
  expect_equal(dim(filtered_vol), dim(vol))
})

test_that("bilateral_filter handles missing mask", {
  vol <- NeuroVol(array(rnorm(20^3), c(20, 20, 20)), NeuroSpace(c(20, 20, 20)))

  filtered_vol <- bilateral_filter(vol, spatial_sigma = 2, intensity_sigma = 1, window = 1)

  expect_s4_class(filtered_vol, "NeuroVol")
  expect_equal(dim(filtered_vol), dim(vol))
  expect_identical(space(filtered_vol), space(vol))
})
