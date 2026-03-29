## Coverage tests for R/neurovol.R
## Targets zero-coverage and low-coverage code paths not exercised by
## test-neurovol.R.

library(neuroim2)

# ---------------------------------------------------------------------------
# Constructors
# ---------------------------------------------------------------------------

test_that("DenseNeuroVol constructor works with 3D array", {
  sp  <- NeuroSpace(c(5L, 5L, 5L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(125), c(5, 5, 5)), sp)
  expect_s4_class(vol, "DenseNeuroVol")
  expect_equal(dim(vol), c(5L, 5L, 5L))
})

test_that("DenseNeuroVol constructor works with 1D vector", {
  sp  <- NeuroSpace(c(5L, 5L, 5L), c(1, 1, 1))
  vol <- DenseNeuroVol(rnorm(125), sp)
  expect_s4_class(vol, "DenseNeuroVol")
  expect_equal(dim(vol), c(5L, 5L, 5L))
})

test_that("DenseNeuroVol constructor uses indices to place sparse data", {
  sp      <- NeuroSpace(c(5L, 5L, 5L), c(1, 1, 1))
  indices <- c(1L, 10L, 50L)
  vol     <- DenseNeuroVol(c(1, 2, 3), sp, indices = indices)
  expect_s4_class(vol, "DenseNeuroVol")
  expect_equal(vol@.Data[indices], c(1, 2, 3))
  expect_equal(vol@.Data[-indices], rep(0, prod(dim(sp)) - 3))
})

test_that("SparseNeuroVol constructor works", {
  sp      <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  indices <- c(1L, 5L, 100L, 500L)
  svol    <- SparseNeuroVol(c(1.1, 2.2, 3.3, 4.4), sp, indices = indices)
  expect_s4_class(svol, "SparseNeuroVol")
  # dense round-trip
  arr <- as.array(svol)
  expect_equal(dim(arr), c(10L, 10L, 10L))
  expect_equal(arr[indices], c(1.1, 2.2, 3.3, 4.4), tolerance = 1e-7)
})

test_that("SparseNeuroVol constructor rejects mismatched data/indices length", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  expect_error(SparseNeuroVol(c(1, 2), sp, indices = c(1L, 2L, 3L)))
})

test_that("SparseNeuroVol constructor rejects out-of-range indices", {
  sp <- NeuroSpace(c(5L, 5L, 5L), c(1, 1, 1))
  expect_error(SparseNeuroVol(1, sp, indices = prod(dim(sp)) + 1L))
})

test_that("LogicalNeuroVol constructor works with logical array", {
  sp  <- NeuroSpace(c(4L, 4L, 4L), c(1, 1, 1))
  arr <- array(c(TRUE, FALSE), c(4, 4, 4))
  lv  <- LogicalNeuroVol(arr, sp)
  expect_s4_class(lv, "LogicalNeuroVol")
  expect_equal(sum(lv), sum(arr))
})

test_that("LogicalNeuroVol constructor coerces numeric array to logical", {
  sp  <- NeuroSpace(c(3L, 3L, 3L), c(1, 1, 1))
  arr <- array(c(0, 1, 2), c(3, 3, 3))   # non-zero treated as TRUE
  lv  <- LogicalNeuroVol(arr, sp)
  expect_s4_class(lv, "LogicalNeuroVol")
  expect_equal(sum(lv), sum(arr > 0))
})

test_that("LogicalNeuroVol constructor works with vector + indices", {
  sp  <- NeuroSpace(c(5L, 5L, 5L), c(1, 1, 1))
  idx <- c(1L, 3L, 7L)
  lv  <- LogicalNeuroVol(rep(TRUE, 3), sp, indices = idx)
  expect_s4_class(lv, "LogicalNeuroVol")
  expect_equal(sum(lv), 3L)
  flat <- as.vector(lv@.Data)
  expect_true(all(flat[idx]))
  expect_true(all(!flat[-idx]))
})

test_that("NeuroVol() factory dispatches to DenseNeuroVol", {
  sp  <- NeuroSpace(c(5L, 5L, 5L), c(1, 1, 1))
  vol <- NeuroVol(array(rnorm(125), c(5, 5, 5)), sp)
  expect_s4_class(vol, "DenseNeuroVol")
})

# ---------------------------------------------------------------------------
# I/O: write_vol / read_vol round-trip
# ---------------------------------------------------------------------------

test_that("write_vol and read_vol round-trip for DenseNeuroVol (.nii)", {
  sp  <- NeuroSpace(c(5L, 5L, 5L), c(2, 2, 2))
  vol <- DenseNeuroVol(array(rnorm(125), c(5, 5, 5)), sp)
  tmp <- tempfile(fileext = ".nii")
  on.exit(unlink(tmp))
  write_vol(vol, tmp)
  vol2 <- read_vol(tmp)
  expect_equal(dim(vol), dim(vol2))
  expect_equal(as.array(vol), as.array(vol2), tolerance = 1e-5)
})

test_that("write_vol and read_vol round-trip for DenseNeuroVol (.nii.gz)", {
  sp  <- NeuroSpace(c(5L, 5L, 5L), c(2, 2, 2))
  vol <- DenseNeuroVol(array(rnorm(125), c(5, 5, 5)), sp)
  tmp <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(tmp))
  write_vol(vol, tmp)
  vol2 <- read_vol(tmp)
  expect_equal(dim(vol), dim(vol2))
  expect_equal(as.array(vol), as.array(vol2), tolerance = 1e-5)
})

# ---------------------------------------------------------------------------
# NeuroVolSource and load_data
# ---------------------------------------------------------------------------

test_that("NeuroVolSource creates valid source from file", {
  f   <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  src <- NeuroVolSource(f, 1)
  expect_s4_class(src, "NeuroVolSource")
  vol <- load_data(src)
  expect_s4_class(vol, "DenseNeuroVol")
  expect_equal(dim(vol), c(64L, 64L, 25L))
})

test_that("read_vol with index exceeding 4D volume count errors", {
  f <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  # file has 4 volumes; index 5 should error
  expect_error(read_vol(f, index = 5))
})

# ---------------------------------------------------------------------------
# Slice extraction
# ---------------------------------------------------------------------------

test_that("slice() extracts 2D slice along axis 3 with correct dims", {
  vol <- make_vol(c(6L, 7L, 8L))
  sl  <- slice(vol, 3L, 3L)
  expect_s4_class(sl, "NeuroSlice")
  expect_equal(dim(sl), c(6L, 7L))
})

test_that("slice() extracts 2D slice along axis 1", {
  vol <- make_vol(c(6L, 7L, 8L))
  sl  <- slice(vol, 2L, 1L)
  expect_s4_class(sl, "NeuroSlice")
  expect_equal(dim(sl), c(7L, 8L))
})

test_that("slice() extracts 2D slice along axis 2", {
  vol <- make_vol(c(6L, 7L, 8L))
  sl  <- slice(vol, 4L, 2L)
  expect_s4_class(sl, "NeuroSlice")
  expect_equal(dim(sl), c(6L, 8L))
})

# ---------------------------------------------------------------------------
# slices() iterator
# ---------------------------------------------------------------------------

test_that("slices() returns deflist of correct length", {
  vol <- make_vol(c(4L, 5L, 6L))
  sls <- slices(vol)
  expect_equal(length(sls), 6L)
  expect_s4_class(sls[[1]], "NeuroSlice")
})

# ---------------------------------------------------------------------------
# map / mapf
# ---------------------------------------------------------------------------

test_that("mapf with Kernel applies smoothing and preserves space", {
  vol  <- make_vol(c(10L, 10L, 10L))
  kern <- Kernel(c(3, 3, 3), vdim = spacing(vol))
  vol2 <- mapf(vol, kern)
  expect_equal(space(vol), space(vol2))
})

# ---------------------------------------------------------------------------
# as.mask
# ---------------------------------------------------------------------------

test_that("as.mask with no indices converts positives to TRUE", {
  vol  <- make_vol(c(4L, 4L, 4L))
  msk  <- as.mask(vol)
  expect_s4_class(msk, "LogicalNeuroVol")
  expect_equal(sum(msk), sum(vol > 0))
})

test_that("as.mask with indices marks only those indices TRUE", {
  sp   <- NeuroSpace(c(4L, 4L, 4L), c(1, 1, 1))
  vol  <- DenseNeuroVol(rep(1, 64), sp)
  idx  <- c(1L, 5L, 10L)
  msk  <- as.mask(vol, indices = idx)
  expect_s4_class(msk, "LogicalNeuroVol")
  expect_equal(sum(msk), 3L)
})

# ---------------------------------------------------------------------------
# as.sparse / as.numeric / as.vector / as.array for SparseNeuroVol
# ---------------------------------------------------------------------------

test_that("as.sparse creates SparseNeuroVol from DenseNeuroVol", {
  vol  <- make_vol()
  idx  <- which(vol > 0)
  svol <- as.sparse(vol, mask = idx)
  expect_s4_class(svol, "SparseNeuroVol")
})

test_that("as.numeric returns dense vector of correct length", {
  sp   <- NeuroSpace(c(5L, 5L, 5L), c(1, 1, 1))
  idx  <- c(1L, 3L, 10L)
  svol <- SparseNeuroVol(c(1.1, 2.2, 3.3), sp, indices = idx)
  v    <- as.numeric(svol)
  expect_equal(length(v), prod(dim(svol)))
})

test_that("as.vector returns dense vector of correct length", {
  sp   <- NeuroSpace(c(5L, 5L, 5L), c(1, 1, 1))
  idx  <- c(1L, 3L, 10L)
  svol <- SparseNeuroVol(c(1.1, 2.2, 3.3), sp, indices = idx)
  v    <- as.vector(svol)
  expect_equal(length(v), prod(dim(svol)))
})

test_that("as.array for SparseNeuroVol has correct dims and values", {
  sp      <- NeuroSpace(c(5L, 5L, 5L), c(1, 1, 1))
  indices <- c(1L, 10L, 100L)
  vals    <- c(7, 8, 9)
  svol    <- SparseNeuroVol(vals, sp, indices = indices)
  arr     <- as.array(svol)
  expect_equal(dim(arr), c(5L, 5L, 5L))
  expect_equal(arr[indices], vals, tolerance = 1e-7)
  expect_equal(sum(arr[-indices]), 0)
})

# ---------------------------------------------------------------------------
# values() methods
# ---------------------------------------------------------------------------

test_that("values() on DenseNeuroVol returns underlying array", {
  vol <- make_vol()
  v   <- values(vol)
  expect_equal(v, vol@.Data)
})

test_that("values() on SparseNeuroVol returns dense numeric vector", {
  sp   <- NeuroSpace(c(5L, 5L, 5L), c(1, 1, 1))
  idx  <- c(1L, 3L, 10L)
  svol <- SparseNeuroVol(c(1.1, 2.2, 3.3), sp, indices = idx)
  v    <- values(svol)
  expect_equal(length(v), prod(dim(svol)))
})

# ---------------------------------------------------------------------------
# map_values
# ---------------------------------------------------------------------------

test_that("map_values with list lookup remaps correctly", {
  sp  <- NeuroSpace(c(3L, 3L, 3L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rep(c(0, 1), 27), c(3, 3, 3)), sp)
  out <- map_values(vol, list("0" = 99, "1" = 42))
  expect_true(all(out[vol == 0] == 99))
  expect_true(all(out[vol == 1] == 42))
})

test_that("map_values with matrix lookup remaps correctly", {
  sp  <- NeuroSpace(c(3L, 3L, 3L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rep(c(0, 1), 27), c(3, 3, 3)), sp)
  lu  <- rbind(c(0, 99), c(1, 42))
  out <- map_values(vol, lu)
  expect_true(all(out[vol == 0] == 99))
  expect_true(all(out[vol == 1] == 42))
})

# ---------------------------------------------------------------------------
# coord_to_grid
# ---------------------------------------------------------------------------

test_that("coord_to_grid dispatches on NeuroVol with matrix and returns 3 values", {
  vol  <- make_vol(c(10L, 10L, 10L))
  crd  <- matrix(c(0, 0, 0), nrow = 1)   # origin
  g    <- coord_to_grid(vol, crd)
  expect_equal(length(g), 3L)
})

test_that("coord_to_grid dispatches on NeuroVol with numeric vector and returns 3 values", {
  vol <- make_vol(c(10L, 10L, 10L))
  crd <- c(0, 0, 0)
  g   <- coord_to_grid(vol, crd)
  expect_equal(length(g), 3L)
})

# ---------------------------------------------------------------------------
# patch_set
# ---------------------------------------------------------------------------

test_that("patch_set returns deflist of correct length (no mask)", {
  vol <- make_vol(c(7L, 7L, 7L))
  ps  <- patch_set(vol, c(3, 3, 3))
  expect_true(length(ps) > 0)
})

# ---------------------------------------------------------------------------
# split_fill
# ---------------------------------------------------------------------------

test_that("split_fill applies function per factor level", {
  vol <- make_vol(c(3L, 3L, 3L))
  fac <- factor(rep(1:3, length.out = prod(dim(vol))))
  out <- split_fill(vol, fac, mean)
  expect_s4_class(out, "DenseNeuroVol")
  expect_equal(dim(out), dim(vol))
})

# ---------------------------------------------------------------------------
# concat (single vol -> NeuroVec)
# ---------------------------------------------------------------------------

test_that("concat on single DenseNeuroVol returns 4D NeuroVec with 1 time point", {
  vol <- make_vol()
  v   <- concat(vol)
  expect_s4_class(v, "NeuroVec")
  expect_equal(dim(v)[4], 1L)
})

# ---------------------------------------------------------------------------
# show methods (smoke test — must not error)
# ---------------------------------------------------------------------------

test_that("show DenseNeuroVol does not error", {
  vol <- make_vol()
  expect_output(show(vol))
})

test_that("show SparseNeuroVol does not error", {
  sp   <- NeuroSpace(c(5L, 5L, 5L), c(1, 1, 1))
  idx  <- c(1L, 3L, 10L)
  svol <- SparseNeuroVol(c(1.1, 2.2, 3.3), sp, indices = idx)
  expect_output(show(svol))
})

test_that("show NeuroVol (generic) does not error", {
  vol <- make_vol()
  expect_output(show(vol))
})
