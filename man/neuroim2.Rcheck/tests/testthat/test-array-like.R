context("array_like generics")

library(neuroim2)

make_vol <- function() {
  sp <- NeuroSpace(c(2,2,2))
  DenseNeuroVol(array(1:8, dim=c(2,2,2)), sp)
}

make_vec <- function() {
  sp <- NeuroSpace(c(2,2,2,2))
  DenseNeuroVec(array(1:16, dim=c(2,2,2,2)), sp)
}

test_that("ArrayLike3D respects drop when only i supplied", {
  vol <- make_vol()
  res <- vol[1, drop=FALSE]
  expect_equal(dim(res), c(1,2,2))
  res2 <- vol[1,, , drop=FALSE]
  expect_equal(dim(res2), c(1,2,2))
  res_drop <- vol[1,,]
  expect_true(is.vector(res_drop) || length(dim(res_drop)) < 3)
})

test_that("ArrayLike3D respects drop when j supplied and k missing", {
  vol <- make_vol()
  res <- vol[,1, , drop=FALSE]
  expect_equal(dim(res), c(2,1,2))
})

test_that("ArrayLike4D numeric missing j path returns correct shape", {
  vec <- make_vec()
  res <- vec[1, , , , drop=FALSE]
  expect_equal(dim(res), c(1,2,2,2))
})

test_that("ArrayLike3D matrix indexing validates bounds", {
  vol <- make_vol()
  bad_idx <- matrix(c(3,1,1), ncol=3)
  expect_error(vol[bad_idx, drop=FALSE])
})

test_that("linear_access works for dense objects", {
  vol <- make_vol()
  expect_equal(linear_access(vol, 1:2), c(1,2))
  vec <- make_vec()
  expect_equal(linear_access(vec, 1:3), as.vector(array(1:16, dim=c(2,2,2,2)))[1:3])
})
