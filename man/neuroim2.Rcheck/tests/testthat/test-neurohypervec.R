context("neurohypervec")

library(neuroim2)

make_hvec <- function() {
  # space: 2x2x2, trials=2, features=3
  sp <- NeuroSpace(c(2,2,2,2,3))
  mask_arr <- array(c(1,0,1,0,1,0,0,1), dim=c(2,2,2)) # 4 active voxels
  mask <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2,2,2)))
  nvox <- sum(mask_arr)
  data <- array(seq_len(3*2*nvox), dim=c(3,2,nvox)) # features x trials x voxels
  list(hvec = NeuroHyperVec(data, sp, mask), sp=sp, mask=mask, data=data)
}

test_that("[ subsetting returns correct dims and values with drop", {
  obj <- make_hvec()
  hvec <- obj$hvec

  # subset all spatial voxels, first trial, first feature
  sub <- hvec[,,,1,1, drop=TRUE]
  expect_equal(dim(sub), c(2,2,2))

  # For an active voxel (1,1,1) corresponding to first voxel in mask
  # data ordering: features x trials x voxels
  expected <- hvec@data[1,1,1]
  expect_equal(sub[1,1,1], expected)

  # Non-drop preserves trial/feature axes
  sub_full <- hvec[1:2,1,1,1:2,1:3, drop=FALSE]
  expect_equal(dim(sub_full), c(2,1,1,2,3))
  expect_equal(sub_full[1,1,1,1,1], hvec@data[1,1,1])
})

test_that("series returns zeros for masked-out voxels and data for active ones", {
  obj <- make_hvec(); hvec <- obj$hvec
  active <- series(hvec, 1,1,1)  # masked in
  expect_equal(dim(active), c(3,2))
  expect_true(all(active == hvec@data[,,1]))

  inactive <- series(hvec, 2,1,2)  # masked out (mask value 0)
  expect_true(all(inactive == 0))
})

test_that("linear_access matches bracket retrieval order", {
  obj <- make_hvec(); hvec <- obj$hvec
  # linear index 1 corresponds to spatial voxel 1, trial 1, feature 1
  lin1 <- linear_access(hvec, 1)
  expect_equal(lin1, hvec@data[1,1,1])

  # choose an active voxel: spatial index 5 (i=1,j=1,k=2) maps to lookup index 3
  spatial_idx <- 5
  trial_idx <- 2
  feature_idx <- 2
  lin_idx <- spatial_idx + (trial_idx - 1) * 8 + (feature_idx - 1) * (8 * 2)
  lin_val <- linear_access(hvec, lin_idx)
  expect_equal(lin_val, hvec@data[feature_idx, trial_idx, 3])
})
