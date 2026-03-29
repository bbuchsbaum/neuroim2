# Regression tests for performance optimizations.
# Ensures fast paths produce identical results to generic paths.

test_that("series(DenseNeuroVec, matrix) matches voxel-by-voxel extraction", {
  vec <- make_vec(c(8, 8, 8), ntime = 10)
  coords <- matrix(c(1,1,1, 4,4,4, 8,8,8), ncol = 3, byrow = TRUE)

  # Fast path (DenseNeuroVec, matrix method)
  fast <- series(vec, coords)

  # Manual extraction via individual voxels (known correct)
  manual <- sapply(1:nrow(coords), function(r) {
    vec[coords[r,1], coords[r,2], coords[r,3], ]
  })

  expect_equal(fast, manual)
})

test_that("series(DenseNeuroVec, matrix) rejects invalid inputs", {
  vec <- make_vec(c(4, 4, 4), ntime = 5)

  # Out-of-bounds
  expect_error(series(vec, matrix(c(5, 1, 1), ncol = 3)))
  expect_error(series(vec, matrix(c(0, 1, 1), ncol = 3)))

  # Fractional
  expect_error(series(vec, matrix(c(1.5, 1, 1), ncol = 3)), "whole-number")

  # Wrong column count
  expect_error(series(vec, matrix(1:4, ncol = 2)), "3 columns")
})

test_that("series(DenseNeuroVec, integer) matches matrix path", {
  vec <- make_vec(c(8, 8, 8), ntime = 10)
  idx <- c(1L, 50L, 100L, 500L)
  coords <- arrayInd(idx, dim(vec)[1:3])

  s_int <- series(vec, idx)
  s_mat <- series(vec, coords)

  expect_equal(s_int, s_mat)
})

test_that("linear_access(DenseNeuroVec) returns correct values for both types", {
  vec <- make_vec(c(4, 4, 4), ntime = 3)
  idx_num <- c(1, 10, 50, 100)
  idx_int <- as.integer(idx_num)

  expect_equal(linear_access(vec, idx_num), vec@.Data[idx_num])
  expect_equal(linear_access(vec, idx_int), vec@.Data[idx_int])
})

test_that("SparseNeuroVec [ extraction is correct", {
  svec <- make_sparse_vec(c(8, 8, 8), ntime = 5, frac = 0.5)
  sub <- svec[2:4, 2:4, 2:4, ]
  expect_equal(length(dim(sub)), 4)
  expect_equal(dim(sub)[1:3], c(3L, 3L, 3L))
})

test_that("SparseNeuroVec as.matrix round-trips correctly", {
  svec <- make_sparse_vec(c(6, 6, 6), ntime = 4, frac = 0.4)
  mat <- as.matrix(svec)
  expect_equal(nrow(mat), prod(dim(svec)[1:3]))
  expect_equal(ncol(mat), dim(svec)[4])
  idx <- indices(svec)
  expect_true(all(mat[-idx, ] == 0))
})
