library(testthat)
library(neuroim2)

test_that("mask works for DenseNeuroVol", {
  # Create a dense volume
  dims <- c(10, 10, 10)
  vol <- NeuroVol(array(rnorm(prod(dims)), dims), NeuroSpace(dims))
  
  m <- mask(vol)
  
  expect_s4_class(m, "LogicalNeuroVol")
  expect_equal(dim(m), dims)
  expect_true(all(m@.Data))  # All values should be TRUE
})

test_that("mask works for DenseNeuroVec", {
  # Create a dense 4D vector
  dims <- c(10, 10, 10, 5)
  data <- array(rnorm(prod(dims)), dims)
  vec <- DenseNeuroVec(data, NeuroSpace(dims))
  
  m <- mask(vec)
  
  expect_s4_class(m, "LogicalNeuroVol")
  expect_equal(dim(m), dims[1:3])  # Mask should be 3D
  expect_true(all(m@.Data))  # All values should be TRUE
})

test_that("mask works for LogicalNeuroVol", {
  # Create a LogicalNeuroVol
  dims <- c(8, 8, 8)
  logical_vol <- LogicalNeuroVol(array(runif(prod(dims)) > 0.5, dims), 
                                 NeuroSpace(dims))
  
  m <- mask(logical_vol)
  
  expect_s4_class(m, "LogicalNeuroVol")
  expect_equal(dim(m), dims)
  expect_true(all(as.array(m)))  # When used as data, mask should be all TRUE
})

test_that("mask works for SparseNeuroVec", {
  # Create a sparse vector with explicit mask
  dims <- c(10, 10, 10)
  mask_array <- array(runif(prod(dims)) > 0.6, dims)  # ~40% sparse
  mask_vol <- LogicalNeuroVol(mask_array, NeuroSpace(dims))
  
  n_time <- 3
  sparse_data <- matrix(rnorm(sum(mask_array) * n_time), 
                       nrow = n_time, 
                       ncol = sum(mask_array))
  
  space_4d <- NeuroSpace(c(dims, n_time))
  svec <- SparseNeuroVec(sparse_data, space_4d, mask_vol)
  
  m <- mask(svec)
  
  expect_s4_class(m, "LogicalNeuroVol")
  expect_equal(dim(m), dims)
  expect_equal(m@.Data, mask_array)  # Should return the stored mask
})

test_that("mask works for ClusteredNeuroVol", {
  # Create a simple clustered volume
  dims <- c(10, 10, 10)
  mask_array <- array(FALSE, dims)
  # Create two small clusters
  mask_array[2:4, 2:4, 2:4] <- TRUE
  mask_array[6:8, 6:8, 6:8] <- TRUE
  
  mask_vol <- LogicalNeuroVol(mask_array, NeuroSpace(dims))
  
  # Create cluster assignments
  n_voxels <- sum(mask_array)
  clusters <- rep(1:2, length.out = n_voxels)
  
  cvol <- ClusteredNeuroVol(mask_vol, clusters)
  
  m <- mask(cvol)
  
  expect_s4_class(m, "LogicalNeuroVol")
  expect_equal(dim(m), dims)
  expect_equal(m@.Data, mask_array)  # Should return the stored mask
})

test_that("mask works for NeuroHyperVec", {
  # Create a NeuroHyperVec
  dims <- c(8, 8, 8)
  mask_array <- array(runif(prod(dims)) > 0.7, dims)  # ~30% sparse
  mask_vol <- LogicalNeuroVol(mask_array, NeuroSpace(dims))
  
  # Create data: trials x features x voxels
  n_trials <- 2
  n_features <- 3
  n_voxels <- sum(mask_array)
  
  # Create data as features x trials x voxels array (NeuroHyperVec expects this order)
  data_array <- array(rnorm(n_features * n_trials * n_voxels),
                     dim = c(n_features, n_trials, n_voxels))
  
  # NeuroHyperVec needs (data, space, mask)
  space_hyper <- NeuroSpace(c(dims, n_trials, n_features))
  hvec <- NeuroHyperVec(data_array, space_hyper, mask_vol)
  
  m <- mask(hvec)
  
  expect_s4_class(m, "LogicalNeuroVol")
  expect_equal(dim(m), dims)
  expect_equal(m@.Data, mask_array)  # Should return the stored mask
})

test_that("mask works for NeuroSlice", {
  # Create a 2D slice
  dims <- c(10, 10)
  slice_data <- array(rnorm(prod(dims)), dims)
  
  # Create appropriate NeuroSpace for 2D slice
  space_2d <- NeuroSpace(dims, spacing = c(1, 1), origin = c(0, 0))
  nslice <- NeuroSlice(slice_data, space_2d)
  
  m <- mask(nslice)
  
  expect_s4_class(m, "LogicalNeuroVol")
  # Mask should be 3D with one dimension = 1
  expect_equal(length(dim(m)), 3)
  expect_true(all(m@.Data))  # All values should be TRUE
})

test_that("mask works for FileBackedNeuroVec", {
  skip_if_not(file.exists(system.file("extdata", "global_mask_v4.nii", 
                                      package = "neuroim2")))
  
  # Use an example file if available
  test_file <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  
  if (file.exists(test_file)) {
    # Try to create a FileBackedNeuroVec
    # This might fail if the file is not 4D, so wrap in try
    fbvec <- try(FileBackedNeuroVec(test_file), silent = TRUE)
    
    if (!inherits(fbvec, "try-error")) {
      m <- mask(fbvec)
      
      expect_s4_class(m, "LogicalNeuroVol")
      expect_equal(dim(m), dim(fbvec)[1:3])
      expect_true(all(m@.Data))  # Should be all TRUE for dense data
    }
  }
})

test_that("mask consistency between sparse and dense representations", {
  # Create a sparse representation
  dims <- c(10, 10, 10)
  mask_array <- array(runif(prod(dims)) > 0.5, dims)
  mask_vol <- LogicalNeuroVol(mask_array, NeuroSpace(dims))
  
  sparse_data <- rnorm(sum(mask_array))
  space_4d <- NeuroSpace(c(dims, 1))
  svec <- SparseNeuroVec(matrix(sparse_data, nrow = 1), space_4d, mask_vol)
  
  # Get mask from sparse representation
  sparse_mask <- mask(svec)
  
  # Create equivalent dense representation
  dense_array <- array(0, dims)
  dense_array[mask_array] <- sparse_data
  dense_vol <- NeuroVol(dense_array, NeuroSpace(dims))
  
  # Dense mask should be all TRUE
  dense_mask <- mask(dense_vol)
  
  expect_s4_class(sparse_mask, "LogicalNeuroVol")
  expect_s4_class(dense_mask, "LogicalNeuroVol")
  
  # Sparse mask should match the original mask
  expect_equal(sparse_mask@.Data, mask_array)
  
  # Dense mask should be all TRUE
  expect_true(all(dense_mask@.Data))
})

test_that("mask preserves spatial information", {
  # Test that spatial metadata is preserved
  dims <- c(12, 14, 16)
  spacing_vals <- c(2, 2.5, 3)
  origin_vals <- c(10, 20, 30)
  
  space <- NeuroSpace(dims, spacing = spacing_vals, origin = origin_vals)
  vol <- NeuroVol(array(rnorm(prod(dims)), dims), space)
  
  m <- mask(vol)
  
  expect_equal(dim(m), dims)
  expect_equal(spacing(m), spacing_vals)
  expect_equal(origin(m), origin_vals)
})

test_that("mask works with BigNeuroVec if available", {
  skip_if_not_installed("bigstatsr")
  
  # Create a BigNeuroVec
  dims <- c(5, 5, 5)
  mask_array <- array(runif(prod(dims)) > 0.5, dims)
  mask_vol <- LogicalNeuroVol(mask_array, NeuroSpace(dims))
  
  n_time <- 2
  data <- matrix(rnorm(sum(mask_array) * n_time), 
                nrow = n_time,
                ncol = sum(mask_array))
  
  space_4d <- NeuroSpace(c(dims, n_time))
  bvec <- BigNeuroVec(data, space_4d, mask_vol)
  
  m <- mask(bvec)
  
  expect_s4_class(m, "LogicalNeuroVol")
  expect_equal(dim(m), dims)
  expect_equal(m@.Data, mask_array)
})