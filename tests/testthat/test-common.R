## Tests for functions in R/common.R
## Targets zero-coverage functions: normalize_mask, data-type lookups,
## grid/index helpers, .getRStorage, .niftiExt, matrixToQuatern/quaternToMatrix,
## and split_reduce for NeuroVec.

library(neuroim2)

# ---------------------------------------------------------------------------
# normalize_mask
# ---------------------------------------------------------------------------

test_that("normalize_mask returns all-TRUE array for NULL", {
  res <- neuroim2:::normalize_mask(NULL, c(3L, 3L, 3L))
  expect_true(is.array(res))
  expect_equal(dim(res), c(3L, 3L, 3L))
  expect_true(all(res))
})

test_that("normalize_mask handles LogicalNeuroVol", {
  sp  <- NeuroSpace(c(4L, 4L, 4L), c(1, 1, 1))
  m   <- array(c(rep(TRUE, 32), rep(FALSE, 32)), c(4L, 4L, 4L))
  lv  <- LogicalNeuroVol(m, sp)
  res <- neuroim2:::normalize_mask(lv, c(4L, 4L, 4L))
  expect_true(is.array(res))
  expect_equal(dim(res), c(4L, 4L, 4L))
  # as.array() may attach dimnames; compare values only
  expect_equal(as.vector(res), as.vector(m))
})

test_that("normalize_mask handles numeric NeuroVol (threshold > 0)", {
  sp  <- NeuroSpace(c(3L, 3L, 3L), c(1, 1, 1))
  dat <- c(rep(0, 13), rep(2, 14))      # 13 zeros, 14 positive
  dv  <- DenseNeuroVol(array(dat, c(3L, 3L, 3L)), sp)
  res <- neuroim2:::normalize_mask(dv, c(3L, 3L, 3L))
  expect_equal(sum(res), 14L)
  expect_true(is.logical(res))
})

test_that("normalize_mask handles integer index vector", {
  D   <- c(5L, 5L, 5L)
  idx <- c(1L, 3L, 10L, 25L)
  res <- neuroim2:::normalize_mask(idx, D)
  expect_equal(dim(res), D)
  expect_true(is.logical(res))
  expect_true(res[1, 1, 1])
  expect_false(res[2, 1, 1])
  expect_equal(sum(res), 4L)
})

test_that("normalize_mask handles flat logical vector of length prod(D)", {
  D   <- c(2L, 3L, 4L)
  lv  <- rep(c(TRUE, FALSE), length.out = prod(D))
  res <- neuroim2:::normalize_mask(lv, D)
  expect_equal(dim(res), D)
  expect_equal(as.vector(res), lv)
})

test_that("normalize_mask errors on LogicalNeuroVol dimension mismatch", {
  sp  <- NeuroSpace(c(4L, 4L, 4L), c(1, 1, 1))
  lv  <- LogicalNeuroVol(array(TRUE, c(4L, 4L, 4L)), sp)
  expect_error(neuroim2:::normalize_mask(lv, c(3L, 3L, 3L)))
})

test_that("normalize_mask errors on NeuroVol dimension mismatch", {
  sp  <- NeuroSpace(c(4L, 4L, 4L), c(1, 1, 1))
  dv  <- DenseNeuroVol(array(1.0, c(4L, 4L, 4L)), sp)
  expect_error(neuroim2:::normalize_mask(dv, c(3L, 3L, 3L)))
})

test_that("normalize_mask errors on unrecognized input", {
  expect_error(neuroim2:::normalize_mask("not_a_mask", c(3L, 3L, 3L)))
})

# ---------------------------------------------------------------------------
# .getDataStorage
# ---------------------------------------------------------------------------

test_that(".getDataStorage maps known NIfTI codes", {
  expect_equal(neuroim2:::.getDataStorage(0),  "UNKNOWN")
  expect_equal(neuroim2:::.getDataStorage(2),  "UBYTE")
  expect_equal(neuroim2:::.getDataStorage(4),  "SHORT")
  expect_equal(neuroim2:::.getDataStorage(8),  "INT")
  expect_equal(neuroim2:::.getDataStorage(16), "FLOAT")
  expect_equal(neuroim2:::.getDataStorage(64), "DOUBLE")
})

test_that(".getDataStorage errors on unknown code", {
  expect_error(neuroim2:::.getDataStorage(999))
})

# ---------------------------------------------------------------------------
# .getDataCode
# ---------------------------------------------------------------------------

test_that(".getDataCode maps known storage names to codes", {
  expect_equal(neuroim2:::.getDataCode("UNKNOWN"), 0L)
  expect_equal(neuroim2:::.getDataCode("FLOAT"),  16L)
  expect_equal(neuroim2:::.getDataCode("DOUBLE"), 64L)
  expect_equal(neuroim2:::.getDataCode("SHORT"),   4L)
})

test_that(".getDataCode errors on unknown type", {
  expect_error(neuroim2:::.getDataCode("COMPLEX"))
})

# ---------------------------------------------------------------------------
# .getDataSize
# ---------------------------------------------------------------------------

test_that(".getDataSize maps known types to byte sizes", {
  expect_equal(neuroim2:::.getDataSize("FLOAT"),   4L)
  expect_equal(neuroim2:::.getDataSize("DOUBLE"),  8L)
  expect_equal(neuroim2:::.getDataSize("SHORT"),   2L)
  expect_equal(neuroim2:::.getDataSize("INT"),     4L)
  expect_equal(neuroim2:::.getDataSize("BINARY"),  1L)
})

test_that(".getDataSize errors on unknown type", {
  expect_error(neuroim2:::.getDataSize("COMPLEX"))
})

# ---------------------------------------------------------------------------
# .getRStorage
# ---------------------------------------------------------------------------

test_that(".getRStorage returns 'integer' for integer-like types", {
  expect_equal(neuroim2:::.getRStorage("SHORT"),   "integer")
  expect_equal(neuroim2:::.getRStorage("INTEGER"), "integer")
  expect_equal(neuroim2:::.getRStorage("BYTE"),    "integer")
  expect_equal(neuroim2:::.getRStorage("UBYTE"),   "integer")
  expect_equal(neuroim2:::.getRStorage("LONG"),    "integer")
  expect_equal(neuroim2:::.getRStorage("BINARY"),  "integer")
})

test_that(".getRStorage returns 'double' for float types", {
  expect_equal(neuroim2:::.getRStorage("FLOAT"),  "double")
  expect_equal(neuroim2:::.getRStorage("DOUBLE"), "double")
})

test_that(".getRStorage is case-insensitive", {
  expect_equal(neuroim2:::.getRStorage("float"),  "double")
  expect_equal(neuroim2:::.getRStorage("short"),  "integer")
})

test_that(".getRStorage errors on unknown type", {
  expect_error(neuroim2:::.getRStorage("GARBAGE"))
})

# ---------------------------------------------------------------------------
# .niftiExt
# ---------------------------------------------------------------------------

test_that(".niftiExt returns correct extensions for nifti-single", {
  ext <- neuroim2:::.niftiExt("nifti-single")
  expect_equal(ext$header, "nii")
  expect_equal(ext$data,   "nii")
})

test_that(".niftiExt returns correct extensions for nifti-pair", {
  ext <- neuroim2:::.niftiExt("nifti-pair")
  expect_equal(ext$header, "hdr")
  expect_equal(ext$data,   "img")
})

test_that(".niftiExt returns correct extensions for nifti-gz", {
  ext <- neuroim2:::.niftiExt("nifti-gz")
  expect_equal(ext$header, "nii.gz")
  expect_equal(ext$data,   "nii.gz")
})

test_that(".niftiExt errors on unknown filetype", {
  expect_error(neuroim2:::.niftiExt("afni"))
  expect_error(neuroim2:::.niftiExt("unknown"))
})

# ---------------------------------------------------------------------------
# .gridToIndex3D and .indexToGrid roundtrip
# ---------------------------------------------------------------------------

test_that(".gridToIndex3D and .indexToGrid roundtrip", {
  set.seed(42)
  dims <- c(5L, 6L, 7L)
  idx  <- sort(sample.int(prod(dims), 20))
  grid <- neuroim2:::.indexToGrid(idx, dims)
  idx2 <- neuroim2:::.gridToIndex3D(dims, grid)
  expect_equal(as.integer(idx2), idx)
})

test_that(".gridToIndex3D accepts a single-row voxel vector", {
  dims <- c(5L, 5L, 5L)
  # voxel (1,1,1) should map to index 1
  idx <- neuroim2:::.gridToIndex3D(dims, c(1L, 1L, 1L))
  expect_equal(as.integer(idx), 1L)
})

test_that(".gridToIndex3D errors when dimensions != 3", {
  expect_error(neuroim2:::.gridToIndex3D(c(5L, 6L), matrix(c(1, 2, 3), 1, 3)))
})

test_that(".gridToIndex3D errors when voxmat has wrong number of columns", {
  expect_error(neuroim2:::.gridToIndex3D(c(5L, 6L, 7L), matrix(1:4, 1, 4)))
})

test_that(".indexToGrid errors on index <= 0", {
  expect_error(neuroim2:::.indexToGrid(0L, c(3L, 3L, 3L)))
})

test_that(".indexToGrid errors on out-of-bounds index", {
  expect_error(neuroim2:::.indexToGrid(28L, c(3L, 3L, 3L)))
})

# ---------------------------------------------------------------------------
# matrixToQuatern / quaternToMatrix roundtrip
# ---------------------------------------------------------------------------

test_that("matrixToQuatern returns quaternion of length 3 and valid qfac", {
  mat <- diag(c(2, 3, 4, 1))
  mat[1:3, 4] <- c(10, 20, 30)
  res <- matrixToQuatern(mat)
  expect_length(res$quaternion, 3)
  expect_true(res$qfac %in% c(-1, 1))
})

test_that("matrixToQuatern and quaternToMatrix roundtrip for identity rotation", {
  # Pure isotropic scaling — rotation part is identity
  mat <- diag(c(2, 3, 4, 1))
  mat[1:3, 4] <- c(10, 20, 30)
  res  <- matrixToQuatern(mat)
  mat2 <- quaternToMatrix(res$quaternion, mat[1:3, 4], c(2, 3, 4), res$qfac)
  expect_equal(mat2[1:3, 1:3], mat[1:3, 1:3], tolerance = 1e-10)
  expect_equal(mat2[1:3, 4],   mat[1:3, 4])
})

test_that("matrixToQuatern and quaternToMatrix roundtrip for arbitrary rotation", {
  # Build a rotation matrix via QR decomposition
  set.seed(7)
  raw   <- matrix(rnorm(9), 3, 3)
  R     <- qr.Q(qr(raw))
  if (det(R) < 0) R[, 3] <- -R[, 3]   # ensure proper rotation
  vox   <- c(1.5, 1.5, 2.0)
  # affine: scale columns by voxel size
  mat   <- diag(4)
  mat[1:3, 1:3] <- R * rep(vox, each = 3)
  mat[1:3, 4]   <- c(5, -3, 2)

  res  <- matrixToQuatern(mat)
  mat2 <- quaternToMatrix(res$quaternion, mat[1:3, 4], vox, res$qfac)
  expect_equal(mat2[1:3, 1:3], mat[1:3, 1:3], tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# split_reduce for NeuroVec — split-by-voxel path
# ---------------------------------------------------------------------------

test_that("split_reduce(NeuroVec, factor, FUN) split-by-voxel returns correct shape", {
  set.seed(1)
  sp  <- NeuroSpace(c(4L, 4L, 4L, 10L), c(1, 1, 1))
  dat <- array(rnorm(4 * 4 * 4 * 10), c(4, 4, 4, 10))
  vec <- DenseNeuroVec(dat, sp)

  n_vox <- prod(c(4L, 4L, 4L))   # 64
  fac   <- factor(rep(1:4, length.out = n_vox))
  out   <- split_reduce(vec, fac, mean)

  expect_equal(nrow(out), 4L)    # 4 factor levels
  expect_equal(ncol(out), 10L)   # 10 time points
  expect_equal(rownames(out), as.character(1:4))
})

test_that("split_reduce(NeuroVec, factor) missing-FUN uses rowMeans", {
  set.seed(2)
  sp  <- NeuroSpace(c(3L, 3L, 3L, 8L), c(1, 1, 1))
  dat <- array(rnorm(3 * 3 * 3 * 8), c(3, 3, 3, 8))
  vec <- DenseNeuroVec(dat, sp)

  n_vox <- prod(c(3L, 3L, 3L))   # 27
  fac   <- factor(rep(c("A", "B", "C"), length.out = n_vox))
  out   <- split_reduce(vec, fac)

  expect_equal(nrow(out), 3L)
  expect_equal(ncol(out), 8L)
  expect_equal(rownames(out), c("A", "B", "C"))
})

test_that("split_reduce(NeuroVec, factor) split-by-time: fac length == n_time dispatches to matrix path", {
  set.seed(3)
  sp  <- NeuroSpace(c(3L, 3L, 3L, 6L), c(1, 1, 1))
  dat <- array(rnorm(3 * 3 * 3 * 6), c(3, 3, 3, 6))
  vec <- DenseNeuroVec(dat, sp)

  # factor over time dimension (length == dim[4] == 6)
  # as.matrix(vec) is voxels x time (27 x 6); the matrix split_reduce
  # method requires fac length == nrow(m) == 27, so fac length 6 is invalid
  fac <- factor(c("a", "a", "b", "b", "c", "c"))
  expect_error(split_reduce(vec, fac, mean))
})

test_that("split_reduce(NeuroVec) errors when fac length is wrong", {
  set.seed(4)
  sp  <- NeuroSpace(c(3L, 3L, 3L, 5L), c(1, 1, 1))
  dat <- array(rnorm(3 * 3 * 3 * 5), c(3, 3, 3, 5))
  vec <- DenseNeuroVec(dat, sp)
  # wrong length: neither n_vox (27) nor n_time (5)
  fac <- factor(rep("x", 10))
  expect_error(split_reduce(vec, fac, mean))
})
