library(testthat)
library(neuroim2)

test_that("downsample with factor works correctly", {
  # Create a test 4D image
  data <- array(rnorm(64 * 64 * 32 * 10), dim = c(64, 64, 32, 10))
  space <- NeuroSpace(dim = c(64, 64, 32, 10), 
                      origin = c(0, 0, 0),
                      spacing = c(2, 2, 2))
  nvec <- DenseNeuroVec(data, space)
  
  # Downsample by factor 0.5
  nvec_down <- downsample(nvec, factor = 0.5)
  
  expect_s4_class(nvec_down, "DenseNeuroVec")
  expect_equal(dim(nvec_down), c(32, 32, 16, 10))
  expect_equal(dim(nvec_down)[4], 10) # Time dimension preserved
  expect_equal(spacing(nvec_down)[1:3], c(4, 4, 4))
})

test_that("downsample with spacing works correctly", {
  data <- array(rnorm(60 * 60 * 30 * 8), dim = c(60, 60, 30, 8))
  space <- NeuroSpace(dim = c(60, 60, 30, 8), 
                      origin = c(0, 0, 0),
                      spacing = c(3, 3, 3))
  nvec <- DenseNeuroVec(data, space)
  
  # Downsample to spacing of 6mm
  nvec_down <- downsample(nvec, spacing = c(6, 6, 6))
  
  expect_s4_class(nvec_down, "DenseNeuroVec")
  expect_equal(dim(nvec_down), c(30, 30, 15, 8))
  expect_equal(dim(nvec_down)[4], 8) # Time dimension preserved
  
  # Check that new spacing is approximately what we requested
  new_spacing <- spacing(nvec_down)[1:3]
  expect_true(all(abs(new_spacing - c(6, 6, 6)) < 0.1))
})

test_that("downsample with outdim works correctly", {
  data <- array(rnorm(64 * 64 * 32 * 5), dim = c(64, 64, 32, 5))
  space <- NeuroSpace(dim = c(64, 64, 32, 5), 
                      origin = c(0, 0, 0),
                      spacing = c(2, 2, 2))
  nvec <- DenseNeuroVec(data, space)
  
  # Downsample to specific dimensions
  nvec_down <- downsample(nvec, outdim = c(32, 32, 16))
  
  expect_s4_class(nvec_down, "DenseNeuroVec")
  expect_equal(dim(nvec_down), c(32, 32, 16, 5))
  expect_equal(dim(nvec_down)[4], 5) # Time dimension preserved
})

test_that("downsample preserves aspect ratio warning", {
  data <- array(rnorm(64 * 64 * 32 * 5), dim = c(64, 64, 32, 5))
  space <- NeuroSpace(dim = c(64, 64, 32, 5), 
                      origin = c(0, 0, 0),
                      spacing = c(2, 2, 2))
  nvec <- DenseNeuroVec(data, space)
  
  # Try to downsample with non-uniform scaling (should warn)
  expect_warning(
    nvec_down <- downsample(nvec, outdim = c(32, 32, 8)),
    "aspect ratio"
  )
})

test_that("downsample with factor = 1 returns same dimensions", {
  data <- array(rnorm(32 * 32 * 16 * 4), dim = c(32, 32, 16, 4))
  space <- NeuroSpace(dim = c(32, 32, 16, 4), 
                      origin = c(0, 0, 0),
                      spacing = c(2, 2, 2))
  nvec <- DenseNeuroVec(data, space)
  
  nvec_down <- downsample(nvec, factor = 1.0)
  
  expect_equal(dim(nvec_down), dim(nvec))
})

test_that("downsample with different factors per dimension", {
  data <- array(rnorm(64 * 64 * 32 * 6), dim = c(64, 64, 32, 6))
  space <- NeuroSpace(dim = c(64, 64, 32, 6), 
                      origin = c(0, 0, 0),
                      spacing = c(2, 2, 2))
  nvec <- DenseNeuroVec(data, space)
  
  # Different factors for each spatial dimension
  nvec_down <- downsample(nvec, factor = c(0.5, 0.5, 0.25))
  
  expect_s4_class(nvec_down, "DenseNeuroVec")
  expect_equal(dim(nvec_down), c(32, 32, 8, 6))
  expect_equal(dim(nvec_down)[4], 6) # Time dimension preserved
})

test_that("downsample errors with invalid parameters", {
  data <- array(rnorm(32 * 32 * 16 * 4), dim = c(32, 32, 16, 4))
  space <- NeuroSpace(dim = c(32, 32, 16, 4), 
                      origin = c(0, 0, 0),
                      spacing = c(2, 2, 2))
  nvec <- DenseNeuroVec(data, space)
  
  # Should error if no parameters specified
  expect_error(downsample(nvec), "Exactly one")
  
  # Should error if multiple parameters specified
  expect_error(downsample(nvec, factor = 0.5, spacing = c(4, 4, 4)), "Exactly one")
  
  # Should error with invalid factor
  expect_error(downsample(nvec, factor = 0), "between 0 .* and 1")
  expect_error(downsample(nvec, factor = 1.5), "between 0 .* and 1")
  
  # Should error with wrong outdim length
  expect_error(downsample(nvec, outdim = c(16, 16)), "exactly 3 values")
})

test_that("downsample handles small images correctly", {
  # Very small image
  data <- array(rnorm(4 * 4 * 4 * 2), dim = c(4, 4, 4, 2))
  space <- NeuroSpace(dim = c(4, 4, 4, 2), 
                      origin = c(0, 0, 0),
                      spacing = c(2, 2, 2))
  nvec <- DenseNeuroVec(data, space)
  
  nvec_down <- downsample(nvec, factor = 0.5)
  
  expect_s4_class(nvec_down, "DenseNeuroVec")
  expect_equal(dim(nvec_down), c(2, 2, 2, 2))
  expect_equal(dim(nvec_down)[4], 2) # Time dimension preserved
})

test_that("downsampled values are reasonable averages", {
  # Create a simple test case with known values
  data <- array(0, dim = c(4, 4, 4, 1))
  # Set a 2x2x2 cube to 8
  data[1:2, 1:2, 1:2, 1] <- 8
  
  space <- NeuroSpace(dim = c(4, 4, 4, 1), 
                      origin = c(0, 0, 0),
                      spacing = c(1, 1, 1))
  nvec <- DenseNeuroVec(data, space)
  
  nvec_down <- downsample(nvec, factor = 0.5)
  
  # The first voxel should be the average of the 2x2x2 cube = 8
  expect_equal(as.array(nvec_down)[1, 1, 1, 1], 8, tolerance = 0.01)
  # Other voxels should be 0
  expect_equal(as.array(nvec_down)[2, 2, 2, 1], 0, tolerance = 0.01)
})

test_that("downsample handles NaN and Inf values", {
  data <- array(1, dim = c(4, 4, 4, 2))
  # Add some NaN and Inf values
  data[1, 1, 1, 1] <- NaN
  data[2, 2, 2, 1] <- Inf
  data[3, 3, 3, 1] <- -Inf
  
  space <- NeuroSpace(dim = c(4, 4, 4, 2), 
                      origin = c(0, 0, 0),
                      spacing = c(1, 1, 1))
  nvec <- DenseNeuroVec(data, space)
  
  # Should not crash
  nvec_down <- downsample(nvec, factor = 0.5)
  
  expect_s4_class(nvec_down, "DenseNeuroVec")
  # Result should have finite values where input had finite values
  result <- as.array(nvec_down)
  expect_true(any(is.finite(result)))
})

test_that("downsample validates spacing values", {
  data <- array(1, dim = c(8, 8, 8, 2))
  space <- NeuroSpace(dim = c(8, 8, 8, 2), 
                      origin = c(0, 0, 0),
                      spacing = c(1, 1, 1))
  nvec <- DenseNeuroVec(data, space)
  
  # Should error with negative spacing
  expect_error(downsample(nvec, spacing = c(-1, 2, 2)), "positive")
  
  # Should error with zero spacing
  expect_error(downsample(nvec, spacing = c(0, 2, 2)), "positive")
  
  # Should error with wrong length spacing
  expect_error(downsample(nvec, spacing = c(2, 2)), "length 3")
})

test_that("downsample validates method parameter", {
  data <- array(1, dim = c(8, 8, 8, 2))
  space <- NeuroSpace(dim = c(8, 8, 8, 2), 
                      origin = c(0, 0, 0),
                      spacing = c(1, 1, 1))
  nvec <- DenseNeuroVec(data, space)
  
  # Should error with unsupported method
  expect_error(downsample(nvec, factor = 0.5, method = "lanczos"), 
               "Only 'box' method")
})

# Tests for 3D DenseNeuroVol downsampling
test_that("downsample 3D vol with factor works correctly", {
  # Create a test 3D volume
  data <- array(rnorm(64 * 64 * 32), dim = c(64, 64, 32))
  space <- NeuroSpace(dim = c(64, 64, 32), 
                      origin = c(0, 0, 0),
                      spacing = c(2, 2, 2))
  vol <- DenseNeuroVol(data, space)
  
  # Downsample by factor 0.5
  vol_down <- downsample(vol, factor = 0.5)
  
  expect_s4_class(vol_down, "DenseNeuroVol")
  expect_equal(dim(vol_down), c(32, 32, 16))
  expect_equal(spacing(vol_down), c(4, 4, 4))
})

test_that("downsample 3D vol with spacing works correctly", {
  data <- array(rnorm(60 * 60 * 30), dim = c(60, 60, 30))
  space <- NeuroSpace(dim = c(60, 60, 30), 
                      origin = c(0, 0, 0),
                      spacing = c(3, 3, 3))
  vol <- DenseNeuroVol(data, space)
  
  # Downsample to spacing of 6mm
  vol_down <- downsample(vol, spacing = c(6, 6, 6))
  
  expect_s4_class(vol_down, "DenseNeuroVol")
  expect_equal(dim(vol_down), c(30, 30, 15))
  
  # Check that new spacing is approximately what we requested
  new_spacing <- spacing(vol_down)
  expect_true(all(abs(new_spacing - c(6, 6, 6)) < 0.1))
})

test_that("downsample 3D vol with outdim works correctly", {
  data <- array(rnorm(64 * 64 * 32), dim = c(64, 64, 32))
  space <- NeuroSpace(dim = c(64, 64, 32), 
                      origin = c(0, 0, 0),
                      spacing = c(2, 2, 2))
  vol <- DenseNeuroVol(data, space)
  
  # Downsample to specific dimensions
  vol_down <- downsample(vol, outdim = c(32, 32, 16))
  
  expect_s4_class(vol_down, "DenseNeuroVol")
  expect_equal(dim(vol_down), c(32, 32, 16))
})

test_that("downsample 3D vol preserves aspect ratio warning", {
  data <- array(rnorm(64 * 64 * 32), dim = c(64, 64, 32))
  space <- NeuroSpace(dim = c(64, 64, 32), 
                      origin = c(0, 0, 0),
                      spacing = c(2, 2, 2))
  vol <- DenseNeuroVol(data, space)
  
  # Try to downsample with non-uniform scaling (should warn)
  expect_warning(
    vol_down <- downsample(vol, outdim = c(32, 32, 8)),
    "aspect ratio"
  )
})

test_that("downsample 3D vol with factor = 1 returns same dimensions", {
  data <- array(rnorm(32 * 32 * 16), dim = c(32, 32, 16))
  space <- NeuroSpace(dim = c(32, 32, 16), 
                      origin = c(0, 0, 0),
                      spacing = c(2, 2, 2))
  vol <- DenseNeuroVol(data, space)
  
  vol_down <- downsample(vol, factor = 1.0)
  
  expect_equal(dim(vol_down), dim(vol))
})

test_that("downsample 3D vol with different factors per dimension", {
  data <- array(rnorm(64 * 64 * 32), dim = c(64, 64, 32))
  space <- NeuroSpace(dim = c(64, 64, 32), 
                      origin = c(0, 0, 0),
                      spacing = c(2, 2, 2))
  vol <- DenseNeuroVol(data, space)
  
  # Different factors for each spatial dimension
  vol_down <- downsample(vol, factor = c(0.5, 0.5, 0.25))
  
  expect_s4_class(vol_down, "DenseNeuroVol")
  expect_equal(dim(vol_down), c(32, 32, 8))
})

test_that("downsample 3D vol errors with invalid parameters", {
  data <- array(rnorm(32 * 32 * 16), dim = c(32, 32, 16))
  space <- NeuroSpace(dim = c(32, 32, 16), 
                      origin = c(0, 0, 0),
                      spacing = c(2, 2, 2))
  vol <- DenseNeuroVol(data, space)
  
  # Should error if no parameters specified
  expect_error(downsample(vol), "Exactly one")
  
  # Should error if multiple parameters specified
  expect_error(downsample(vol, factor = 0.5, spacing = c(4, 4, 4)), "Exactly one")
  
  # Should error with invalid factor
  expect_error(downsample(vol, factor = 0), "between 0 .* and 1")
  expect_error(downsample(vol, factor = 1.5), "between 0 .* and 1")
  
  # Should error with wrong outdim length
  expect_error(downsample(vol, outdim = c(16, 16)), "exactly 3 values")
})

test_that("downsample 3D vol handles small volumes correctly", {
  # Very small volume
  data <- array(rnorm(4 * 4 * 4), dim = c(4, 4, 4))
  space <- NeuroSpace(dim = c(4, 4, 4), 
                      origin = c(0, 0, 0),
                      spacing = c(2, 2, 2))
  vol <- DenseNeuroVol(data, space)
  
  vol_down <- downsample(vol, factor = 0.5)
  
  expect_s4_class(vol_down, "DenseNeuroVol")
  expect_equal(dim(vol_down), c(2, 2, 2))
})

test_that("3D downsampled values are reasonable averages", {
  # Create a simple test case with known values
  data <- array(0, dim = c(4, 4, 4))
  # Set a 2x2x2 cube to 8
  data[1:2, 1:2, 1:2] <- 8
  
  space <- NeuroSpace(dim = c(4, 4, 4), 
                      origin = c(0, 0, 0),
                      spacing = c(1, 1, 1))
  vol <- DenseNeuroVol(data, space)
  
  vol_down <- downsample(vol, factor = 0.5)
  
  # The first voxel should be the average of the 2x2x2 cube = 8
  expect_equal(as.array(vol_down)[1, 1, 1], 8, tolerance = 0.01)
  # Other voxels should be 0
  expect_equal(as.array(vol_down)[2, 2, 2], 0, tolerance = 0.01)
})

test_that("downsample 3D vol handles NaN and Inf values", {
  data <- array(1, dim = c(4, 4, 4))
  # Add some NaN and Inf values
  data[1, 1, 1] <- NaN
  data[2, 2, 2] <- Inf
  data[3, 3, 3] <- -Inf
  
  space <- NeuroSpace(dim = c(4, 4, 4), 
                      origin = c(0, 0, 0),
                      spacing = c(1, 1, 1))
  vol <- DenseNeuroVol(data, space)
  
  # Should not crash
  vol_down <- downsample(vol, factor = 0.5)
  
  expect_s4_class(vol_down, "DenseNeuroVol")
  # Result should have finite values where input had finite values
  result <- as.array(vol_down)
  expect_true(any(is.finite(result)))
})

test_that("downsample 3D vol validates spacing values", {
  data <- array(1, dim = c(8, 8, 8))
  space <- NeuroSpace(dim = c(8, 8, 8), 
                      origin = c(0, 0, 0),
                      spacing = c(1, 1, 1))
  vol <- DenseNeuroVol(data, space)
  
  # Should error with negative spacing
  expect_error(downsample(vol, spacing = c(-1, 2, 2)), "positive")
  
  # Should error with zero spacing
  expect_error(downsample(vol, spacing = c(0, 2, 2)), "positive")
  
  # Should error with wrong length spacing
  expect_error(downsample(vol, spacing = c(2, 2)), "length 3")
})

test_that("downsample 3D vol validates method parameter", {
  data <- array(1, dim = c(8, 8, 8))
  space <- NeuroSpace(dim = c(8, 8, 8), 
                      origin = c(0, 0, 0),
                      spacing = c(1, 1, 1))
  vol <- DenseNeuroVol(data, space)
  
  # Should error with unsupported method
  expect_error(downsample(vol, factor = 0.5, method = "lanczos"), 
               "Only 'box' method")
})
