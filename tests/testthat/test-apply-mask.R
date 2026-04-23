library(testthat)
library(neuroim2)

test_that("apply_mask zeroes voxels outside the mask for DenseNeuroVol", {
  sp <- NeuroSpace(c(4, 4, 4), spacing = c(2, 2, 2), origin = c(10, 20, 30))
  arr <- array(seq_len(prod(dim(sp))), dim(sp))
  vol <- NeuroVol(arr, sp)

  mask_arr <- array(FALSE, dim(sp))
  mask_arr[c(1, 5, 9)] <- TRUE
  mask_vol <- LogicalNeuroVol(mask_arr, sp)

  out <- apply_mask(vol, mask_vol)

  expect_s4_class(out, "DenseNeuroVol")
  expect_equal(dim(out), dim(vol))
  expect_equal(spacing(out), spacing(vol))
  expect_equal(origin(out), origin(vol))
  expect_equal(as.array(out)[mask_arr], arr[mask_arr])
  expect_true(all(as.array(out)[!mask_arr] == 0))
})

test_that("apply_mask preserves logical type for LogicalNeuroVol", {
  sp <- NeuroSpace(c(3, 3, 3))
  x_arr <- array(FALSE, dim(sp))
  x_arr[c(1, 2, 3)] <- TRUE
  x <- LogicalNeuroVol(x_arr, sp)

  mask_arr <- array(FALSE, dim(sp))
  mask_arr[c(2, 3, 4)] <- TRUE

  out <- apply_mask(x, mask_arr)

  expect_s4_class(out, "LogicalNeuroVol")
  expect_equal(as.vector(as.array(out)), as.vector(x_arr & mask_arr))
})

test_that("apply_mask zeroes voxels outside the mask for DenseNeuroVec", {
  dims <- c(3, 3, 3, 2)
  sp <- NeuroSpace(dims, spacing = c(1, 1, 1), origin = c(0, 0, 0))
  arr <- array(seq_len(prod(dims)), dims)
  vec <- DenseNeuroVec(arr, sp, label = "run", volume_labels = c("a", "b"))

  mask_arr <- array(FALSE, dims[1:3])
  mask_arr[c(1, 3, 5)] <- TRUE

  out <- apply_mask(vec, mask_arr)
  out_arr <- as.array(out)
  expected <- arr
  for (tt in seq_len(dims[4])) {
    vol <- expected[, , , tt]
    vol[!mask_arr] <- 0
    expected[, , , tt] <- vol
  }

  expect_s4_class(out, "DenseNeuroVec")
  expect_equal(volume_labels(out), c("a", "b"))
  expect_equal(as.vector(out_arr), as.vector(expected))
})

test_that("apply_mask intersects masks for sparse vectors", {
  dims <- c(4, 4, 4, 2)
  sp <- NeuroSpace(dims)

  src_mask <- array(FALSE, dims[1:3])
  src_mask[c(1, 2, 3, 4)] <- TRUE
  data <- matrix(seq_len(sum(src_mask) * dims[4]), nrow = dims[4], ncol = sum(src_mask))
  svec <- SparseNeuroVec(data, sp, src_mask, label = "sparse")

  keep_mask <- array(FALSE, dims[1:3])
  keep_mask[c(1, 4)] <- TRUE

  out <- apply_mask(svec, keep_mask)
  out_arr <- as.array(out)
  expected_mask <- src_mask & keep_mask
  expected <- as.array(svec)
  for (tt in seq_len(dims[4])) {
    vol <- expected[, , , tt]
    vol[!expected_mask] <- 0
    expected[, , , tt] <- vol
  }

  expect_s4_class(out, "SparseNeuroVec")
  expect_equal(as.vector(as.array(mask(out))), as.vector(expected_mask))
  expect_equal(as.vector(out_arr), as.vector(expected))
})

test_that("apply_mask validates mask geometry for NeuroVol masks", {
  sp <- NeuroSpace(c(3, 3, 3), spacing = c(2, 2, 2))
  vol <- NeuroVol(array(1, dim(sp)), sp)
  bad_mask <- NeuroVol(array(1, dim(sp)), NeuroSpace(c(3, 3, 3), spacing = c(1, 1, 1)))

  expect_error(apply_mask(vol, bad_mask), "Spatial geometry")
})
