# First, make sure to load the necessary packages
library(testthat)

# Load or source the code containing the NeuroHyperVec and NeuroSpace classes
# Assume the classes are defined in 'neuroimaging_classes.R'
# source('neuroimaging_classes.R')

# Start writing the test suite
context("NeuroHyperVec Class Tests")

test_that("NeuroHyperVec constructor works correctly with valid inputs", {
  # Define dimensions
  spatial_dims <- c(10L, 10L, 10L)
  num_trials <- 5L
  num_features <- 3L
  
  # Create a NeuroSpace object
  space <- NeuroSpace(
    dim = c(spatial_dims, num_trials, num_features),
    spacing = c(2.0, 2.0, 2.0),
    origin = c(-10.0, -10.0, -10.0)
  )
  
  # Create a mask with 20% of voxels active
  set.seed(123)
  mask_data <- array(runif(prod(spatial_dims)) < 0.2, dim = spatial_dims)
  mask <- LogicalNeuroVol(mask_data, NeuroSpace(spatial_dims))
  
  # Generate random data for active voxels
  num_voxels <- sum(mask_data)
  data_array <- array(rnorm(num_features * num_trials * num_voxels),
                      dim = c(num_features, num_trials, num_voxels))
  
  # Create NeuroHyperVec object
  hvec <- NeuroHyperVec(data = data_array, space = space, mask = mask)
  
  # Check that the object is created correctly
  expect_s4_class(hvec, "NeuroHyperVec")
  expect_equal(dim(hvec@data), c(num_features, num_trials, num_voxels))
  expect_equal(hvec@space, space)
  expect_equal(hvec@mask, mask)
  expect_length(hvec@lookup_map, prod(spatial_dims))
  expect_equal(sum(hvec@lookup_map > 0), num_voxels)
})

test_that("NeuroHyperVec constructor throws errors with invalid inputs", {
  # Define dimensions
  spatial_dims <- c(10L, 10L, 10L)
  num_trials <- 5L
  num_features <- 3L
  
  # Create a NeuroSpace object
  space <- NeuroSpace(
    dim = c(spatial_dims, num_trials, num_features),
    spacing = c(2.0, 2.0, 2.0),
    origin = c(-10.0, -10.0, -10.0)
  )
  
  # Create a mask
  mask_data <- array(TRUE, dim = spatial_dims)
  mask <- LogicalNeuroVol(mask_data, NeuroSpace(spatial_dims))
  
  # Generate data with incorrect dimensions
  data_array <- array(rnorm(100), dim = c(2, 2, 2))  # Wrong dimensions
  
  # Attempt to create NeuroHyperVec object with incorrect data dimensions
  expect_error(
    NeuroHyperVec(data = data_array, space = space, mask = mask),
    "Data array dimensions \\[2 x 2 x 2\\] do not match expected \\[3 x 5 x 1000\\]"
  )
  
  # Attempt to create NeuroHyperVec object with invalid mask
  invalid_mask <- array(TRUE, dim = c(5, 5))  # Wrong dimensions
  expect_error(
    NeuroHyperVec(data = data_array, space = space, mask = invalid_mask),
    "'mask' must be a 3D logical array"
  )
})

test_that("series method retrieves correct data", {
  # Define dimensions
  spatial_dims <- c(5L, 5L, 5L)
  num_trials <- 2L
  num_features <- 2L
  
  # Create a NeuroSpace object
  space <- NeuroSpace(
    dim = c(spatial_dims, num_trials, num_features),
    spacing = c(1.0, 1.0, 1.0),
    origin = c(0.0, 0.0, 0.0)
  )
  
  # Create a mask with all voxels active
  mask_data <- array(TRUE, dim = spatial_dims)
  mask <- LogicalNeuroVol(mask_data, NeuroSpace(spatial_dims))
  
  # Generate sequential data for testing
  num_voxels <- sum(mask_data)
  data_array <- array(1:(num_features * num_trials * num_voxels),
                      dim = c(num_features, num_trials, num_voxels))
  
  # Create NeuroHyperVec object
  hvec <- NeuroHyperVec(data = data_array, space = space, mask = mask)
  
  # Test series method at a specific voxel
  i <- 2L; j <- 2L; k <- 2L
  result <- series(hvec, i, j, k)
  
  # Expected data
  voxel_index <- ((k - 1) * spatial_dims[1] * spatial_dims[2]) +
                 ((j - 1) * spatial_dims[1]) + i
  lookup_index <- hvec@lookup_map[voxel_index]
  expected_data <- hvec@data[,,lookup_index]
  
  expect_equal(result, expected_data)
})

test_that("series method returns zeros for voxels outside mask", {
  # Define dimensions
  spatial_dims <- c(5L, 5L, 5L)
  num_trials <- 2L
  num_features <- 2L
  
  # Create a NeuroSpace object
  space <- NeuroSpace(
    dim = c(spatial_dims, num_trials, num_features),
    spacing = c(1.0, 1.0, 1.0),
    origin = c(0.0, 0.0, 0.0)
  )
  
  # Create a mask with some voxels inactive
  mask_data <- array(FALSE, dim = spatial_dims)
  mask_data[2,2,2] <- TRUE
  mask <- LogicalNeuroVol(mask_data, NeuroSpace(spatial_dims))
  
  # Generate data
  num_voxels <- sum(mask_data)
  data_array <- array(rnorm(num_features * num_trials * num_voxels),
                      dim = c(num_features, num_trials, num_voxels))
  
  # Create NeuroHyperVec object
  hvec <- NeuroHyperVec(data = data_array, space = space, mask = mask)
  
  # Test series method at an inactive voxel
  i <- 3L; j <- 3L; k <- 3L
  result <- series(hvec, i, j, k)
  
  # Expected data is zeros
  expected_data <- array(0, dim = c(num_features, num_trials))
  
  expect_equal(result, expected_data)
})

test_that("linear_access method retrieves correct data", {
  # Define dimensions
  spatial_dims <- c(4L, 4L, 4L)
  num_trials <- 2L
  num_features <- 2L
  
  # Create a NeuroSpace object
  space <- NeuroSpace(
    dim = c(spatial_dims, num_trials, num_features),
    spacing = c(1.0, 1.0, 1.0),
    origin = c(0.0, 0.0, 0.0)
  )
  
  # Create a mask with all voxels active
  mask_data <- array(TRUE, dim = spatial_dims)
  mask <- LogicalNeuroVol(mask_data, NeuroSpace(spatial_dims))
  
  # Generate data
  num_voxels <- sum(mask_data)
  data_array <- array(1:(num_features * num_trials * num_voxels),
                      dim = c(num_features, num_trials, num_voxels))
  
  # Create NeuroHyperVec object
  hvec <- NeuroHyperVec(data = data_array, space = space, mask = mask)
  
  # Total number of elements
  total_elements <- prod(c(spatial_dims, num_trials, num_features))
  
  # Test linear_access at random indices
  set.seed(456)
  indices <- sample(1:total_elements, 10)
  result <- linear_access(hvec, indices)
  
  # Expected data
  # Map indices to (feature_idx, trial_idx, spatial_idx)
  tmp <- indices - 1
  spatial_nels <- prod(spatial_dims)
  spatial_idx <- (tmp %% spatial_nels) + 1
  tmp <- tmp %/% spatial_nels
  trial_idx <- (tmp %% num_trials) + 1
  feature_idx <- (tmp %/% num_trials) + 1
  
  lookup_indices <- hvec@lookup_map[spatial_idx]
  expected_data <- numeric(length(indices))
  expected_data <- hvec@data[cbind(feature_idx, trial_idx, lookup_indices)]
  
  expect_equal(result, expected_data)
})

test_that("linear_access returns zeros for indices outside mask", {
  # Define dimensions
  spatial_dims <- c(4L, 4L, 4L)
  num_trials <- 2L
  num_features <- 2L
  
  # Create a NeuroSpace object
  space <- NeuroSpace(
    dim = c(spatial_dims, num_trials, num_features),
    spacing = c(1.0, 1.0, 1.0),
    origin = c(0.0, 0.0, 0.0)
  )
  
  # Create a mask with some voxels inactive
  mask_data <- array(FALSE, dim = spatial_dims)
  mask_data[1,1,1] <- TRUE
  mask <- LogicalNeuroVol(mask_data, NeuroSpace(spatial_dims))
  
  # Generate data
  num_voxels <- sum(mask_data)
  data_array <- array(rnorm(num_features * num_trials * num_voxels),
                      dim = c(num_features, num_trials, num_voxels))
  
  # Create NeuroHyperVec object
  hvec <- NeuroHyperVec(data = data_array, space = space, mask = mask)
  
  # Indices corresponding to inactive voxels
  total_elements <- prod(c(spatial_dims, num_trials, num_features))
  indices <- c(1, total_elements)  # First and last indices
  
  # Test linear_access
  result <- linear_access(hvec, indices)
  
  # Expected data is zero for the last index (inactive voxel)
  expect_true(result[1] != 0)
  expect_equal(result[2], 0)
})

test_that("Extracting subsets with [ works correctly", {
  # Define dimensions
  spatial_dims <- c(5L, 5L, 5L)
  num_trials <- 2L
  num_features <- 2L
  
  # Create a NeuroSpace object
  space <- NeuroSpace(
    dim = c(spatial_dims, num_trials, num_features),
    spacing = c(1.0, 1.0, 1.0),
    origin = c(0.0, 0.0, 0.0)
  )
  
  # Create a mask with all voxels active
  mask_data <- array(TRUE, dim = spatial_dims)
  mask <- LogicalNeuroVol(mask_data, NeuroSpace(spatial_dims))
  
  # Generate data
  num_voxels <- sum(mask_data)
  data_array <- array(1:(num_features * num_trials * num_voxels),
                      dim = c(num_features, num_trials, num_voxels))
  
  # Create NeuroHyperVec object
  hvec <- NeuroHyperVec(data = data_array, space = space, mask = mask)
  
  # Extract a subset
  subset_data <- hvec[1:2, 1:2, 1:2, 1, 1, drop = TRUE]
  
  # Expected data dimensions
  expect_equal(dim(subset_data), c(2, 2, 2))
})

test_that("show method works without errors", {
  # Define dimensions
  spatial_dims <- c(10L, 10L, 10L)
  num_trials <- 3L
  num_features <- 2L
  
  # Create a NeuroSpace object
  space <- NeuroSpace(
    dim = c(spatial_dims, num_trials, num_features),
    spacing = c(2.0, 2.0, 2.0),
    origin = c(-10.0, -10.0, -10.0)
  )
  
  # Create a mask with 50% of voxels active
  set.seed(789)
  mask_data <- array(runif(prod(spatial_dims)) < 0.5, dim = spatial_dims)
  mask <- LogicalNeuroVol(mask_data, NeuroSpace(spatial_dims))
  
  # Generate data
  num_voxels <- sum(mask_data)
  data_array <- array(rnorm(num_features * num_trials * num_voxels),
                      dim = c(num_features, num_trials, num_voxels))
  
  # Create NeuroHyperVec object
  hvec <- NeuroHyperVec(data = data_array, space = space, mask = mask)
  
  # Capture output of show method
  expect_output(show(hvec))
})

test_that("NeuroSpace is used correctly within NeuroHyperVec", {
  # Define dimensions
  spatial_dims <- c(8L, 8L, 8L)
  num_trials <- 4L
  num_features <- 2L
  
  # Create a NeuroSpace object with extra dimensions
  space <- NeuroSpace(
    dim = c(spatial_dims, num_trials, num_features),
    spacing = c(1.5, 1.5, 1.5),
    origin = c(-12.0, -12.0, -12.0)
  )
  
  # Check that spacing and origin have length 3
  expect_equal(length(space@spacing), 3)
  expect_equal(length(space@origin), 3)
  
  # Create a mask
  mask_data <- array(runif(prod(spatial_dims)) < 0.3, dim = spatial_dims)
  mask <- LogicalNeuroVol(mask_data, NeuroSpace(spatial_dims))
  
  # Generate data
  num_voxels <- sum(mask_data)
  data_array <- array(rnorm(num_features * num_trials * num_voxels),
                      dim = c(num_features, num_trials, num_voxels))
  
  # Create NeuroHyperVec object
  hvec <- NeuroHyperVec(data = data_array, space = space, mask = mask)
  
  # Check that NeuroSpace properties are consistent
  expect_equal(dim(hvec@space), c(spatial_dims, num_trials, num_features))
  expect_equal(hvec@space@spacing, c(1.5, 1.5, 1.5))
  expect_equal(hvec@space@origin, c(-12.0, -12.0, -12.0))
})

test_that("Error handling for invalid indices in series method", {
  # Define dimensions
  spatial_dims <- c(5L, 5L, 5L)
  num_trials <- 2L
  num_features <- 2L
  
  # Create a NeuroSpace object
  space <- NeuroSpace(
    dim = c(spatial_dims, num_trials, num_features),
    spacing = c(1.0, 1.0, 1.0),
    origin = c(0.0, 0.0, 0.0)
  )
  
  # Create a mask
  mask_data <- array(TRUE, dim = spatial_dims)
  mask <- LogicalNeuroVol(mask_data, NeuroSpace(spatial_dims))
  
  # Generate data
  num_voxels <- sum(mask_data)
  data_array <- array(rnorm(num_features * num_trials * num_voxels),
                      dim = c(num_features, num_trials, num_voxels))
  
  # Create NeuroHyperVec object
  hvec <- NeuroHyperVec(data = data_array, space = space, mask = mask)
  
  # Attempt to access invalid indices
  expect_error(
    series(hvec, 6, 1, 1),
    "Indices out of bounds"
  )
  expect_error(
    series(hvec, 1, 6, 1),
    "Indices out of bounds"
  )
  expect_error(
    series(hvec, 1, 1, 6),
    "Indices out of bounds"
  )
})

test_that("Error handling for invalid indices in linear_access method", {
  # Define dimensions
  spatial_dims <- c(4L, 4L, 4L)
  num_trials <- 2L
  num_features <- 2L
  
  # Create a NeuroSpace object
  space <- NeuroSpace(
    dim = c(spatial_dims, num_trials, num_features),
    spacing = c(1.0, 1.0, 1.0),
    origin = c(0.0, 0.0, 0.0)
  )
  
  # Create a mask
  mask_data <- array(TRUE, dim = spatial_dims)
  mask <- LogicalNeuroVol(mask_data, NeuroSpace(spatial_dims))
  
  # Generate data
  num_voxels <- sum(mask_data)
  data_array <- array(rnorm(num_features * num_trials * num_voxels),
                      dim = c(num_features, num_trials, num_voxels))
  
  # Create NeuroHyperVec object
  hvec <- NeuroHyperVec(data = data_array, space = space, mask = mask)
  
  # Total number of elements
  total_elements <- prod(c(spatial_dims, num_trials, num_features))
  
  # Attempt to access invalid indices
  expect_error(
    linear_access(hvec, c(0, total_elements + 1)),
    "indices must be within range of data dimensions"
  )
})