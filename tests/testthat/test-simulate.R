library(testthat)
library(neuroim2)

test_that("simulate_fmri creates correct dimensions", {
  # Create a small test mask
  dims <- c(16, 16, 10)
  mask_array <- array(TRUE, dims)
  mask <- NeuroVol(mask_array, NeuroSpace(dims, c(2, 2, 2)))
  
  n_time <- 50
  sim_data <- simulate_fmri(mask, n_time = n_time, seed = 123)
  
  expect_s4_class(sim_data, "NeuroVec")
  expect_equal(dim(sim_data), c(dims, n_time))
})

test_that("simulate_fmri works with binary mask", {
  # Create a spherical mask
  dims <- c(20, 20, 12)
  mask_array <- array(FALSE, dims)
  center <- dims / 2
  radius_sq <- 36  # 6^2
  
  for (i in 1:dims[1]) {
    for (j in 1:dims[2]) {
      for (k in 1:dims[3]) {
        dist_sq <- sum((c(i, j, k) - center)^2)
        if (dist_sq <= radius_sq) {
          mask_array[i, j, k] <- TRUE
        }
      }
    }
  }
  
  mask <- NeuroVol(mask_array, NeuroSpace(dims, c(3, 3, 3)))
  n_time <- 30
  
  sim_data <- simulate_fmri(mask, n_time = n_time, seed = 456)
  
  expect_s4_class(sim_data, "NeuroVec")
  expect_equal(dim(sim_data), c(dims, n_time))
  
  # Check that values outside mask are zero
  sim_array <- as.array(sim_data)
  for (t in 1:n_time) {
    expect_true(all(sim_array[,,,t][!mask_array] == 0))
  }
})

test_that("simulate_fmri produces temporal autocorrelation", {
  dims <- c(10, 10, 8)
  mask_array <- array(TRUE, dims)
  mask <- NeuroVol(mask_array, NeuroSpace(dims))
  
  n_time <- 100
  sim_data <- simulate_fmri(
    mask, 
    n_time = n_time, 
    ar_mean = 0.6,  # Higher AR for stronger autocorrelation
    ar_sd = 0.05,
    seed = 789
  )
  
  # Check a few random voxels for positive lag-1 autocorrelation
  sim_array <- as.array(sim_data)
  n_check <- 5
  mask_idx <- which(mask_array)
  check_idx <- sample(mask_idx, n_check)
  
  for (idx in check_idx) {
    coords <- arrayInd(idx, dims)
    ts <- sim_array[coords[1], coords[2], coords[3], ]
    
    # Lag-1 autocorrelation should be positive
    if (sd(ts) > 0) {
      ac <- acf(ts, lag.max = 1, plot = FALSE)$acf[2]
      expect_true(ac > 0.2)  # Should have some autocorrelation
    }
  }
})

test_that("simulate_fmri centering works correctly", {
  dims <- c(8, 8, 6)
  mask_array <- array(TRUE, dims)
  mask <- NeuroVol(mask_array, NeuroSpace(dims))
  
  n_time <- 40
  
  # Test with centering
  sim_centered <- simulate_fmri(
    mask, 
    n_time = n_time, 
    return_centered = TRUE,
    seed = 111
  )
  
  sim_array <- as.array(sim_centered)
  # Check a few voxels for near-zero mean
  n_check <- 10
  mask_idx <- which(mask_array)
  check_idx <- sample(mask_idx, n_check)
  
  for (idx in check_idx) {
    coords <- arrayInd(idx, dims)
    ts <- sim_array[coords[1], coords[2], coords[3], ]
    expect_true(abs(mean(ts)) < 0.1)  # Should be close to zero
  }
  
  # Test without centering
  sim_uncentered <- simulate_fmri(
    mask, 
    n_time = n_time, 
    return_centered = FALSE,
    seed = 222
  )
  
  # At least some voxels should have non-zero mean
  sim_array2 <- as.array(sim_uncentered)
  means <- numeric(n_check)
  for (i in 1:n_check) {
    coords <- arrayInd(check_idx[i], dims)
    ts <- sim_array2[coords[1], coords[2], coords[3], ]
    means[i] <- mean(ts)
  }
  expect_true(sd(means) > 0)  # Should have variation in means
})

test_that("simulate_fmri works without optional components", {
  dims <- c(10, 10, 6)
  mask_array <- array(TRUE, dims)
  mask <- NeuroVol(mask_array, NeuroSpace(dims))
  
  n_time <- 20
  
  # Disable global signal and factors
  sim_data <- simulate_fmri(
    mask,
    n_time = n_time,
    global_amp = 0,      # No global signal
    n_factors = 0,       # No latent factors
    seed = 333
  )
  
  expect_s4_class(sim_data, "NeuroVec")
  expect_equal(dim(sim_data), c(dims, n_time))
})

test_that("simulate_fmri reproducible with seed", {
  dims <- c(8, 8, 4)
  mask_array <- array(TRUE, dims)
  mask <- NeuroVol(mask_array, NeuroSpace(dims))
  
  n_time <- 10
  
  sim1 <- simulate_fmri(mask, n_time = n_time, seed = 999)
  sim2 <- simulate_fmri(mask, n_time = n_time, seed = 999)
  
  expect_equal(as.array(sim1), as.array(sim2))
})

test_that("simulate_fmri handles empty mask gracefully", {
  dims <- c(10, 10, 6)
  mask_array <- array(FALSE, dims)  # Empty mask
  mask <- NeuroVol(mask_array, NeuroSpace(dims))
  
  expect_error(
    simulate_fmri(mask, n_time = 10),
    "Mask is empty"
  )
})

test_that("simulate_fmri TR attribute is set", {
  dims <- c(8, 8, 4)
  mask_array <- array(TRUE, dims)
  mask <- NeuroVol(mask_array, NeuroSpace(dims))
  
  TR_val <- 1.5
  sim_data <- simulate_fmri(mask, n_time = 10, TR = TR_val, seed = 123)
  
  expect_equal(attr(sim_data, "TR"), TR_val)
})
