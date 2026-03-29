library(testthat)
library(neuroim2)

test_that("downsample with factor works correctly", {
  skip_on_cran()
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
  skip_on_cran()
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

test_that("downsample updates 4D affine consistently with new spacing", {
  data <- array(rnorm(20 * 30 * 40 * 2), dim = c(20, 30, 40, 2))
  affine <- matrix(c(
    0, 2, 0, -20,
    -2, 0, 0, 10,
    0, 0, 2, 30,
    0, 0, 0, 1
  ), nrow = 4, byrow = TRUE)
  nvec <- DenseNeuroVec(
    data,
    NeuroSpace(dim = c(20, 30, 40, 2), spacing = c(2, 2, 2), trans = affine)
  )

  nvec_down <- downsample(nvec, spacing = c(4, 4, 4))
  expected_trans <- rescale_affine(affine, c(20, 30, 40), c(4, 4, 4), c(10, 15, 20))

  expect_equal(trans(nvec_down), expected_trans, tolerance = 1e-6)
  expect_equal(voxel_sizes(trans(nvec_down)), c(4, 4, 4), tolerance = 1e-6)
})

test_that("downsample with outdim works correctly", {
  skip_on_cran()
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
  skip_on_cran()
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
  skip_on_cran()
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
  skip_on_cran()
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
  skip_on_cran()
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

test_that("downsample updates 3D affine consistently with new spacing", {
  data <- array(rnorm(20 * 30 * 40), dim = c(20, 30, 40))
  affine <- matrix(c(
    0, 0, 1.5, -40,
    0, -1.5, 0, 25,
    1.5, 0, 0, 5,
    0, 0, 0, 1
  ), nrow = 4, byrow = TRUE)
  vol <- DenseNeuroVol(
    data,
    NeuroSpace(dim = c(20, 30, 40), spacing = c(1.5, 1.5, 1.5), trans = affine)
  )

  vol_down <- downsample(vol, spacing = c(3, 3, 3))
  expected_trans <- rescale_affine(affine, c(20, 30, 40), c(3, 3, 3), c(10, 15, 20))

  expect_equal(trans(vol_down), expected_trans, tolerance = 1e-6)
  expect_equal(voxel_sizes(trans(vol_down)), c(3, 3, 3), tolerance = 1e-6)
})

test_that("downsample 3D vol with outdim works correctly", {
  skip_on_cran()
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
  skip_on_cran()
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
  skip_on_cran()
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

make_sparse_neurovec <- function(coords, values, dims = c(4, 4, 4), spacing = c(1, 1, 1),
                                 trans = NULL, label = "") {
  stopifnot(ncol(coords) == 3)
  stopifnot(ncol(values) == nrow(coords))

  full_dims <- c(dims, nrow(values))
  if (is.null(trans)) {
    sp <- NeuroSpace(dim = full_dims, spacing = spacing)
  } else {
    sp <- NeuroSpace(dim = full_dims, spacing = spacing, trans = trans)
  }

  mask_arr <- array(FALSE, dim = dims)
  linear_idx <- 1L +
    (coords[, 1] - 1L) +
    (coords[, 2] - 1L) * dims[1] +
    (coords[, 3] - 1L) * dims[1] * dims[2]
  mask_arr[linear_idx] <- TRUE

  SparseNeuroVec(values, sp, mask_arr, label = label)
}

downsample_sparse_reference <- function(x, new_dims) {
  old_dims <- dim(x)[1:3]
  input_coords <- arrayInd(indices(x), .dim = old_dims)
  output_coords <- cbind(
    floor((input_coords[, 1] - 1) * new_dims[1] / old_dims[1]) + 1L,
    floor((input_coords[, 2] - 1) * new_dims[2] / old_dims[2]) + 1L,
    floor((input_coords[, 3] - 1) * new_dims[3] / old_dims[3]) + 1L
  )

  output_keys <- do.call(paste, c(as.data.frame(output_coords), sep = ":"))
  groups <- split(seq_len(nrow(output_coords)), output_keys)
  first_members <- vapply(groups, `[`, integer(1), 1L)
  kept_coords <- output_coords[first_members, , drop = FALSE]
  linear_idx <- 1L +
    (kept_coords[, 1] - 1L) +
    (kept_coords[, 2] - 1L) * new_dims[1] +
    (kept_coords[, 3] - 1L) * new_dims[1] * new_dims[2]
  order_idx <- order(linear_idx)
  groups <- groups[order_idx]

  out_data <- vapply(groups, function(idx) {
    rowMeans(x@data[, idx, drop = FALSE])
  }, numeric(dim(x)[4]))

  if (!is.matrix(out_data)) {
    out_data <- matrix(out_data, nrow = dim(x)[4], ncol = 1L)
  }

  out_mask <- array(FALSE, dim = new_dims)
  out_mask[linear_idx[order_idx]] <- TRUE

  list(data = out_data, mask = out_mask)
}

test_that("downsample SparseNeuroVec with factor = 1 preserves values, mask, and label", {
  coords <- rbind(
    c(1, 1, 1),
    c(2, 3, 1),
    c(4, 4, 4)
  )
  values <- rbind(
    c(1, 2, 3),
    c(4, 5, 6)
  )
  svec <- make_sparse_neurovec(coords, values, label = "sparse-id")

  out <- downsample(svec, factor = 1)

  expect_s4_class(out, "SparseNeuroVec")
  expect_equal(dim(out), dim(svec))
  expect_equal(as.array(mask(out)), as.array(mask(svec)))
  expect_equal(unname(out@data), unname(svec@data))
  expect_equal(out@label, "sparse-id")
})

test_that("downsample SparseNeuroVec uses missing-aware box means", {
  coords <- rbind(
    c(1, 1, 1),
    c(2, 2, 2),
    c(4, 4, 4)
  )
  values <- rbind(
    c(2, 4, 8),
    c(6, 10, 14)
  )
  svec <- make_sparse_neurovec(coords, values)

  out <- downsample(svec, factor = 0.5)

  expect_equal(dim(out), c(2, 2, 2, 2))
  expect_equal(sum(as.array(mask(out))), 2)
  expect_true(as.array(mask(out))[1, 1, 1])
  expect_true(as.array(mask(out))[2, 2, 2])
  expect_equal(unname(out@data[, 1]), c(3, 8))
  expect_equal(unname(out@data[, 2]), c(8, 14))
})

test_that("downsample SparseNeuroVec differs from dense zero-filled averaging", {
  dims <- c(4, 4, 4, 1)
  data <- array(0, dim = dims)
  data[1, 1, 1, 1] <- 4
  data[2, 2, 2, 1] <- 8

  dense <- DenseNeuroVec(data, NeuroSpace(dim = dims, spacing = c(1, 1, 1)))
  sparse <- make_sparse_neurovec(
    coords = rbind(c(1, 1, 1), c(2, 2, 2)),
    values = matrix(c(4, 8), nrow = 1)
  )

  dense_out <- downsample(dense, factor = 0.5)
  sparse_out <- downsample(sparse, factor = 0.5)

  expect_equal(as.array(dense_out)[1, 1, 1, 1], 1.5)
  expect_equal(sparse_out@data[1, 1], 6)
  expect_false(isTRUE(all.equal(as.array(dense_out)[1, 1, 1, 1], sparse_out@data[1, 1])))
})

test_that("downsample SparseNeuroVec matches grouped reference implementation", {
  coords <- rbind(
    c(1, 1, 1),
    c(1, 2, 2),
    c(2, 2, 2),
    c(3, 3, 4),
    c(4, 4, 3)
  )
  values <- rbind(
    c(2, 4, 6, 8, 10),
    c(1, 3, 5, 7, 9),
    c(0, 2, 4, 6, 8)
  )
  svec <- make_sparse_neurovec(coords, values, dims = c(4, 4, 4), spacing = c(2, 2, 2))

  out <- downsample(svec, factor = 0.5)
  ref <- downsample_sparse_reference(svec, c(2, 2, 2))

  expect_equal(unname(out@data), unname(ref$data))
  expect_equal(dim(as.array(mask(out))), dim(ref$mask))
  expect_equal(as.vector(as.array(mask(out))), as.vector(ref$mask))
})

test_that("downsample SparseNeuroVec updates affine and spacing consistently", {
  affine <- matrix(c(
    0, 2, 0, -20,
    -2, 0, 0, 10,
    0, 0, 2, 30,
    0, 0, 0, 1
  ), nrow = 4, byrow = TRUE)
  coords <- rbind(
    c(1, 1, 1),
    c(2, 2, 2),
    c(3, 3, 3),
    c(4, 4, 4)
  )
  values <- rbind(
    c(1, 2, 3, 4),
    c(5, 6, 7, 8)
  )
  svec <- make_sparse_neurovec(coords, values, spacing = c(2, 2, 2), trans = affine)

  out <- downsample(svec, spacing = c(4, 4, 4))
  expected_trans <- rescale_affine(affine, c(4, 4, 4), c(4, 4, 4), c(2, 2, 2))

  expect_equal(spacing(out), c(4, 4, 4))
  expect_equal(trans(out), expected_trans, tolerance = 1e-6)
})

test_that("downsample SparseNeuroVec accepts equivalent factor, spacing, and outdim targets", {
  coords <- rbind(
    c(1, 1, 1),
    c(2, 2, 2),
    c(3, 3, 3),
    c(4, 4, 4)
  )
  values <- rbind(
    c(1, 2, 3, 4),
    c(5, 6, 7, 8)
  )
  svec <- make_sparse_neurovec(coords, values, spacing = c(2, 2, 2))

  by_factor <- downsample(svec, factor = 0.5)
  by_spacing <- downsample(svec, spacing = c(4, 4, 4))
  by_outdim <- downsample(svec, outdim = c(2, 2, 2))

  expect_equal(by_factor@data, by_spacing@data)
  expect_equal(by_factor@data, by_outdim@data)
  expect_equal(as.array(mask(by_factor)), as.array(mask(by_spacing)))
  expect_equal(as.array(mask(by_factor)), as.array(mask(by_outdim)))
})

test_that("downsample SparseNeuroVec validates method parameter", {
  coords <- rbind(c(1, 1, 1), c(2, 2, 2))
  values <- matrix(c(1, 2), nrow = 1)
  svec <- make_sparse_neurovec(coords, values)

  expect_error(downsample(svec, factor = 0.5, method = "lanczos"),
               "Only 'box' method")
})
