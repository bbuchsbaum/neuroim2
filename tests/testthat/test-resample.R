library(testthat)
library(neuroim2)



# Test the resample(NeuroVol, NeuroVol) function
test_that("resample(NeuroVol, NeuroVol) works correctly", {
  source <- NeuroVol(array(rnorm(64 * 64 * 64), c(64, 64, 64)), NeuroSpace(c(64, 64, 64)))
  target <- NeuroVol(array(rnorm(64 * 64 * 64), c(64, 64, 64)), NeuroSpace(c(64, 64, 64)))
  interpolation <- 3L

  resampled_vol <- resample(source, target, interpolation)

  expect_s4_class(resampled_vol, "NeuroVol")
  expect_equal(dim(resampled_vol), dim(target))
  expect_equal(space(resampled_vol), space(target))
})

# Test the resample(NeuroVol, NeuroSpace) function
test_that("resample(NeuroVol, NeuroSpace) works correctly", {
  source <- NeuroVol(array(rnorm(64 * 64 * 64), c(64, 64, 64)), NeuroSpace(c(64, 64, 64)))
  target_space <- NeuroSpace(c(64, 64, 64))
  interpolation <- 3L

  resampled_vol <- resample(source, target_space, interpolation)

  expect_s3_class(resampled_vol, "NeuroVol")
  expect_equal(dim(resampled_vol), dim(target_space))
  expect_equal(space(resampled_vol), target_space)
})
