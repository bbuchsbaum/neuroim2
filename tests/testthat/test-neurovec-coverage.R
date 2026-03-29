## Coverage tests for R/neurovec.R
## Targets zero-coverage and low-coverage code paths not exercised by
## test-neurovec.R.

library(neuroim2)

# ---------------------------------------------------------------------------
# NeuroVec constructor dispatching
# ---------------------------------------------------------------------------

test_that("NeuroVec() creates DenseNeuroVec from 4D array", {
  sp  <- NeuroSpace(c(5L, 5L, 5L, 10L), c(1, 1, 1))
  arr <- array(rnorm(5 * 5 * 5 * 10), c(5, 5, 5, 10))
  vec <- NeuroVec(arr, sp)
  expect_s4_class(vec, "DenseNeuroVec")
  expect_equal(dim(vec), c(5L, 5L, 5L, 10L))
})

test_that("NeuroVec() with mask creates SparseNeuroVec", {
  sp   <- NeuroSpace(c(5L, 5L, 5L, 10L), c(1, 1, 1))
  arr  <- array(rnorm(5 * 5 * 5 * 10), c(5, 5, 5, 10))
  mask <- array(FALSE, c(5, 5, 5))
  mask[1:10] <- TRUE
  vec  <- NeuroVec(arr, sp, mask = mask)
  expect_s4_class(vec, "SparseNeuroVec")
  expect_equal(dim(vec), c(5L, 5L, 5L, 10L))
})

test_that("NeuroVec() from list of NeuroVol builds correct 4D object", {
  sp3  <- NeuroSpace(c(5L, 5L, 5L), c(1, 1, 1))
  vols <- lapply(seq_len(3), function(i)
    DenseNeuroVol(array(rnorm(125), c(5, 5, 5)), sp3))
  vec  <- NeuroVec(vols)
  expect_s4_class(vec, "DenseNeuroVec")
  expect_equal(dim(vec)[1:3], c(5L, 5L, 5L))
  expect_equal(dim(vec)[4], 3L)
})

# ---------------------------------------------------------------------------
# vec_from_vols
# ---------------------------------------------------------------------------

test_that("vec_from_vols builds NeuroVec from list of vols", {
  sp3  <- NeuroSpace(c(4L, 4L, 4L), c(1, 1, 1))
  vols <- lapply(seq_len(5), function(i)
    DenseNeuroVol(rnorm(64), sp3))
  vec  <- vec_from_vols(vols)
  expect_s4_class(vec, "NeuroVec")
  expect_equal(dim(vec)[4], 5L)
})

test_that("vec_from_vols rejects non-list input", {
  sp3 <- NeuroSpace(c(4L, 4L, 4L), c(1, 1, 1))
  vol <- DenseNeuroVol(rnorm(64), sp3)
  expect_error(vec_from_vols(vol))
})

test_that("vec_from_vols rejects empty list", {
  expect_error(vec_from_vols(list()))
})

# ---------------------------------------------------------------------------
# DenseNeuroVec constructor edge cases
# ---------------------------------------------------------------------------

test_that("DenseNeuroVec constructor handles matrix (voxels x time)", {
  sp  <- NeuroSpace(c(5L, 5L, 5L, 4L), c(1, 1, 1))
  mat <- matrix(rnorm(500), nrow = 125, ncol = 4)
  vec <- DenseNeuroVec(mat, sp)
  expect_s4_class(vec, "DenseNeuroVec")
  expect_equal(dim(vec), c(5L, 5L, 5L, 4L))
})

test_that("DenseNeuroVec constructor handles matrix (time x voxels) with auto-transpose", {
  sp  <- NeuroSpace(c(5L, 5L, 5L, 4L), c(1, 1, 1))
  mat <- matrix(rnorm(500), nrow = 4, ncol = 125)
  vec <- DenseNeuroVec(mat, sp)
  expect_s4_class(vec, "DenseNeuroVec")
  expect_equal(dim(vec), c(5L, 5L, 5L, 4L))
})

# ---------------------------------------------------------------------------
# write_vec / read_vec round-trip
# ---------------------------------------------------------------------------

test_that("write_vec and read_vec round-trip for DenseNeuroVec (.nii)", {
  vec <- make_vec(c(5L, 5L, 5L), ntime = 4L)
  tmp <- tempfile(fileext = ".nii")
  on.exit(unlink(tmp))
  write_vec(vec, tmp)
  vec2 <- read_vec(tmp)
  expect_equal(dim(vec), dim(vec2))
  # compare data only (labels differ between original and read-back)
  expect_equal(vec@.Data, vec2@.Data, tolerance = 1e-5)
})

# ---------------------------------------------------------------------------
# NeuroVecSource
# ---------------------------------------------------------------------------

test_that("NeuroVecSource builds source from 4D file", {
  f   <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  src <- neuroim2:::NeuroVecSource(f)
  expect_s4_class(src, "NeuroVecSource")
})

test_that("NeuroVecSource with indices subset loads correct number of volumes", {
  f   <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  src <- neuroim2:::NeuroVecSource(f, indices = 1:2)
  vec <- load_data(src)
  expect_equal(dim(vec)[4], 2L)
})

# ---------------------------------------------------------------------------
# sub_vector
# ---------------------------------------------------------------------------

test_that("sub_vector extracts time subset correctly", {
  vec <- make_vec(c(5L, 5L, 5L), ntime = 10L)
  sv  <- sub_vector(vec, 1:3)
  expect_s4_class(sv, "NeuroVec")
  expect_equal(dim(sv)[4], 3L)
  expect_equal(dim(sv)[1:3], dim(vec)[1:3])
})

test_that("sub_vector with single index returns 4D with 1 time point", {
  vec <- make_vec(c(4L, 4L, 4L), ntime = 8L)
  sv  <- sub_vector(vec, 5L)
  expect_equal(dim(sv)[4], 1L)
})

# ---------------------------------------------------------------------------
# series with various index types
# ---------------------------------------------------------------------------

test_that("series with integer linear index returns correct time series", {
  vec  <- make_vec(c(5L, 5L, 5L), ntime = 6L)
  idx  <- 1L
  ts   <- series(vec, idx)
  # should equal vec[1,1,1,]
  expect_equal(ts, as.numeric(vec[1, 1, 1, ]))
})

test_that("series with matrix coordinates returns correct result", {
  bv  <- make_vec(c(6L, 6L, 6L), ntime = 4L)
  mat <- rbind(c(1, 1, 1), c(2, 2, 2), c(3, 3, 3))
  r1  <- apply(mat, 1, function(i) series(bv, i[1], i[2], i[3]))
  r2  <- series(bv, mat)
  expect_equal(r1, r2)
})

test_that("series with i,j,k numeric coordinates returns single time series", {
  vec <- make_vec(c(5L, 5L, 5L), ntime = 6L)
  ts  <- series(vec, 2, 3, 4)
  expect_equal(length(ts), 6L)
  expect_equal(ts, as.numeric(vec[2, 3, 4, ]))
})

test_that("series with LogicalNeuroVol mask returns matrix", {
  vec  <- make_vec(c(5L, 5L, 5L), ntime = 4L)
  mask <- make_mask(c(5L, 5L, 5L), frac = 0.4)
  result <- series(vec, mask)
  nmask  <- sum(mask)
  expect_equal(ncol(result), nmask)
  expect_equal(nrow(result), 4L)
})

test_that("series with NeuroVol (non-logical) dispatches via as.logical", {
  vec <- make_vec(c(5L, 5L, 5L), ntime = 4L)
  sp3 <- NeuroSpace(c(5L, 5L, 5L), c(1, 1, 1))
  nvol <- DenseNeuroVol(array(c(0, 1), c(5, 5, 5)), sp3)
  result <- series(vec, nvol)
  expect_equal(ncol(result), sum(nvol > 0))
})

# ---------------------------------------------------------------------------
# series_roi
# ---------------------------------------------------------------------------

test_that("series_roi with matrix returns ROIVec", {
  vec <- make_vec(c(6L, 6L, 6L), ntime = 4L)
  mat <- rbind(c(1, 1, 1), c(2, 2, 2))
  roi <- series_roi(vec, mat)
  expect_s4_class(roi, "ROIVec")
  expect_equal(nrow(coords(roi)), 2L)
})

# ---------------------------------------------------------------------------
# as.matrix / as.list
# ---------------------------------------------------------------------------

test_that("as.matrix returns voxels x time matrix", {
  vec <- make_vec(c(4L, 4L, 4L), ntime = 5L)
  m   <- as.matrix(vec)
  expect_equal(nrow(m), prod(dim(vec)[1:3]))
  expect_equal(ncol(m), dim(vec)[4])
})

test_that("as.list returns list of NeuroVol objects", {
  vec <- make_vec(c(4L, 4L, 4L), ntime = 3L)
  lst <- as.list(vec)
  expect_equal(length(lst), 3L)
  expect_s4_class(lst[[1]], "NeuroVol")
  expect_equal(dim(lst[[1]]), c(4L, 4L, 4L))
})

# ---------------------------------------------------------------------------
# length method
# ---------------------------------------------------------------------------

test_that("length() on NeuroVec returns number of time points", {
  vec <- make_vec(c(5L, 5L, 5L), ntime = 7L)
  expect_equal(length(vec), 7L)
})

# ---------------------------------------------------------------------------
# drop method
# ---------------------------------------------------------------------------

test_that("drop() on NeuroVec with 1 time point returns NeuroVol", {
  sp  <- NeuroSpace(c(5L, 5L, 5L, 1L), c(1, 1, 1))
  arr <- array(rnorm(125), c(5, 5, 5, 1))
  vec <- DenseNeuroVec(arr, sp)
  result <- drop(vec)
  expect_true(inherits(result, "NeuroVol"))
  expect_equal(dim(result), c(5L, 5L, 5L))
})

test_that("drop() on NeuroVec with > 1 time points returns same NeuroVec", {
  vec <- make_vec(c(5L, 5L, 5L), ntime = 3L)
  result <- drop(vec)
  expect_true(inherits(result, "NeuroVec"))
  expect_equal(dim(result)[4], 3L)
})

# ---------------------------------------------------------------------------
# concat NeuroVol + NeuroVec and NeuroVec + NeuroVol
# ---------------------------------------------------------------------------

test_that("concat(NeuroVec, NeuroVol) increases time dimension by 1", {
  vec <- make_vec(c(4L, 4L, 4L), ntime = 3L)
  vol <- make_vol(c(4L, 4L, 4L))
  out <- concat(vec, vol)
  expect_equal(dim(out)[4], 4L)
})

test_that("concat(NeuroVol, NeuroVec) increases time dimension by 1", {
  vol <- make_vol(c(4L, 4L, 4L))
  vec <- make_vec(c(4L, 4L, 4L), ntime = 3L)
  out <- concat(vol, vec)
  expect_equal(dim(out)[4], 4L)
})

# ---------------------------------------------------------------------------
# as.sparse (DenseNeuroVec -> SparseNeuroVec)
# ---------------------------------------------------------------------------

test_that("as.sparse with LogicalNeuroVol mask converts DenseNeuroVec", {
  vec  <- make_vec(c(5L, 5L, 5L), ntime = 4L)
  mask <- make_mask(c(5L, 5L, 5L), frac = 0.5)
  svec <- as.sparse(vec, mask)
  expect_s4_class(svec, "SparseNeuroVec")
  expect_equal(dim(svec)[1:3], c(5L, 5L, 5L))
  expect_equal(dim(svec)[4], 4L)
})

test_that("as.sparse with numeric indices converts DenseNeuroVec", {
  vec <- make_vec(c(5L, 5L, 5L), ntime = 4L)
  idx <- c(1L, 5L, 25L, 50L)
  svec <- as.sparse(vec, idx)
  expect_s4_class(svec, "SparseNeuroVec")
})

# ---------------------------------------------------------------------------
# read_vol_list
# ---------------------------------------------------------------------------

test_that("read_vol_list loads multiple 3D volume files into a NeuroVec", {
  f    <- system.file("extdata", "global_mask2.nii.gz", package = "neuroim2")
  vec  <- read_vol_list(list(f, f))
  expect_s4_class(vec, "NeuroVec")
  expect_equal(dim(vec)[4], 2L)
})

# ---------------------------------------------------------------------------
# show methods (smoke tests)
# ---------------------------------------------------------------------------

test_that("show DenseNeuroVec does not error", {
  vec <- make_vec()
  expect_output(show(vec))
})

test_that("show NeuroVecSource does not error", {
  f   <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  src <- neuroim2:::NeuroVecSource(f)
  expect_output(show(src))
})
