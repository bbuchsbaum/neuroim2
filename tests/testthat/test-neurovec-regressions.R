context("neurovec regressions")

library(neuroim2)
library(RNifti)

test_that("DenseNeuroVec preserves array input shape", {
  arr <- array(1:16, dim = c(2,2,2,2))
  sp  <- NeuroSpace(c(2,2,2,2))

  vec <- DenseNeuroVec(arr, sp)

  expect_equal(dim(vec), dim(arr))
  expect_identical(as(vec, "array"), arr)
})

test_that("series_roi accepts j/k with drop flag", {
  arr <- array(seq_len(24), dim = c(2,2,2,3))
  sp  <- NeuroSpace(c(2,2,2,3))
  vec <- DenseNeuroVec(arr, sp)

  v_drop    <- series_roi(vec, i = 1, j = 1, k = 1)
  v_nodrop  <- series_roi(vec, i = 1, j = 1, k = 1, drop = FALSE)

  expect_equal(v_drop, arr[1,1,1,])
  expect_equal(as.numeric(v_nodrop), arr[1,1,1,])
  expect_true(length(v_nodrop) == 3)
})

test_that("read_vol_list builds correct NeuroVec from files", {
  td <- tempdir()
  sp <- NeuroSpace(c(2,2,2))
  arr1 <- array(1:8, dim = c(2,2,2))
  arr2 <- array(seq(10, 80, by = 10), dim = c(2,2,2))

  f1 <- file.path(td, "test_vol1.nii")
  f2 <- file.path(td, "test_vol2.nii")
  RNifti::writeNifti(arr1, f1)
  RNifti::writeNifti(arr2, f2)

  vec <- read_vol_list(c(f1, f2))
  expect_s4_class(vec, "NeuroVec")
  expect_equal(dim(vec), c(2,2,2,2))

  v1 <- vec[[1]]
  v2 <- vec[[2]]
  expect_equal(v1@.Data, arr1)
  expect_equal(v2@.Data, arr2)
})
