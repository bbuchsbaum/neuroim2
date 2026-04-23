library(testthat)
library(neuroim2)

make_mask_backend_fixture <- function() {
  dims <- c(7, 6, 5, 4)
  sp <- NeuroSpace(dims)
  arr <- array(rnorm(prod(dims), sd = 0.15), dims)

  for (tt in seq_len(dims[4])) {
    arr[2:6, 2:5, 2:4, tt] <- arr[2:6, 2:5, 2:4, tt] + 4 + tt
  }

  dense <- DenseNeuroVec(arr, sp)
  path <- tempfile(fileext = ".nii")
  write_vec(dense, path)

  list(
    path = path,
    dense = dense,
    mmap = read_vec(path, mode = "mmap"),
    filebacked = read_vec(path, mode = "filebacked")
  )
}

test_that("clip_level returns a positive scalar for bright foreground volumes", {
  dims <- c(11, 11, 11)
  sp <- NeuroSpace(dims)
  arr <- array(0, dims)
  arr[3:9, 3:9, 3:9] <- array(seq(60, 110, length.out = 7^3), c(7, 7, 7))
  vol <- NeuroVol(arr, sp)

  lev <- clip_level(vol)

  expect_true(is.numeric(lev))
  expect_length(lev, 1)
  expect_gt(lev, 0)
  expect_lt(lev, max(arr))
})

test_that("clip_level respects integer-valued histograms", {
  dims <- c(11, 11, 11)
  sp <- NeuroSpace(dims)
  arr <- array(as.integer(0), dims)
  arr[3:9, 3:9, 3:9] <- as.integer(array(seq(60, 110, length.out = 7^3), c(7, 7, 7)))
  vol <- NeuroVol(arr, sp)

  lev <- clip_level(vol)

  expect_equal(lev, round(lev), tolerance = 0)
  expect_gt(lev, 0)
})

test_that("clip_level can return a gradual clip map", {
  dims <- c(11, 11, 11)
  sp <- NeuroSpace(dims)
  arr <- array(abs(rnorm(prod(dims))), dims)
  arr[3:9, 3:9, 3:9] <- arr[3:9, 3:9, 3:9] + 20
  vol <- NeuroVol(arr, sp)

  grad <- clip_level(vol, gradual = TRUE)

  expect_s4_class(grad, "DenseNeuroVol")
  expect_equal(dim(grad), dims)
  expect_true(all(as.array(grad) > 0))
})

test_that("clip_level on DenseNeuroVec matches an explicit mean-absolute reference volume", {
  dims <- c(7, 6, 5, 4)
  sp <- NeuroSpace(dims)
  arr <- array(rnorm(prod(dims)), dims)
  arr[2:6, 2:5, 2:4, ] <- arr[2:6, 2:5, 2:4, ] + 8

  vec <- DenseNeuroVec(arr, sp)
  ref <- rowMeans(abs(matrix(arr, nrow = prod(dims[1:3]), ncol = dims[4])))
  ref_vol <- NeuroVol(array(ref, dim = dims[1:3]), drop_dim(sp))

  expect_equal(
    clip_level(vec, representative = "mean_abs", gradual = FALSE),
    clip_level(ref_vol, gradual = FALSE)
  )

  expect_equal(
    as.array(clip_level(vec, representative = "mean_abs", gradual = TRUE)),
    as.array(clip_level(ref_vol, gradual = TRUE))
  )
})

test_that("clip_level on DenseNeuroVec matches an explicit voxelwise median volume", {
  dims <- c(7, 6, 5, 4)
  sp <- NeuroSpace(dims)
  arr <- array(rnorm(prod(dims)), dims)
  arr[2:6, 2:5, 2:4, ] <- arr[2:6, 2:5, 2:4, ] + 8

  vec <- DenseNeuroVec(arr, sp)
  ref <- apply(matrix(arr, nrow = prod(dims[1:3]), ncol = dims[4]), 1L, stats::median)
  ref_vol <- NeuroVol(array(ref, dim = dims[1:3]), drop_dim(sp))

  expect_equal(
    clip_level(vec, representative = "median", gradual = FALSE),
    clip_level(ref_vol, gradual = FALSE)
  )

  expect_equal(
    as.array(clip_level(vec, representative = "median", gradual = TRUE)),
    as.array(clip_level(ref_vol, gradual = TRUE))
  )
})

test_that("mean-absolute representative is invariant to sign flips", {
  dims <- c(7, 6, 5, 4)
  sp <- NeuroSpace(dims)
  arr <- array(rnorm(prod(dims)), dims)
  arr[2:6, 2:5, 2:4, ] <- arr[2:6, 2:5, 2:4, ] + 6

  vec <- DenseNeuroVec(arr, sp)
  flipped <- DenseNeuroVec(-arr, sp)

  expect_equal(
    clip_level(vec, representative = "mean_abs", gradual = FALSE),
    clip_level(flipped, representative = "mean_abs", gradual = FALSE)
  )

  expect_equal(
    as.array(clip_level(vec, representative = "mean_abs", gradual = TRUE)),
    as.array(clip_level(flipped, representative = "mean_abs", gradual = TRUE))
  )

  expect_equal(
    as.array(automask(vec, representative = "mean_abs", gradual = FALSE, peels = 0)),
    as.array(automask(flipped, representative = "mean_abs", gradual = FALSE, peels = 0))
  )
})

test_that("automask returns a logical mask and keeps the largest component", {
  dims <- c(11, 11, 11)
  sp <- NeuroSpace(dims)
  arr <- array(0, dims)
  arr[4:10, 4:10, 4:10] <- array(seq(80, 120, length.out = 7^3), c(7, 7, 7))
  arr[1:2, 1:2, 1:2] <- array(seq(70, 77, length.out = 8), c(2, 2, 2))
  vol <- NeuroVol(arr, sp)

  msk <- automask(vol, gradual = FALSE, peels = 0)
  marr <- as.array(msk)

  expect_s4_class(msk, "LogicalNeuroVol")
  expect_equal(dim(msk), dims)
  expect_true(all(marr[5:9, 5:9, 5:9]))
  expect_false(any(marr[1:2, 1:2, 1:2]))
})

test_that("automask on DenseNeuroVec uses a 3D representative image", {
  dims <- c(11, 11, 11, 3)
  sp <- NeuroSpace(dims)
  arr <- array(0, dims)
  for (tt in seq_len(dims[4])) {
    arr[4:10, 4:10, 4:10, tt] <- array(seq(70 + 5 * tt, 110 + 5 * tt, length.out = 7^3), c(7, 7, 7))
  }

  vec <- DenseNeuroVec(arr, sp)
  msk <- automask(vec, gradual = FALSE, peels = 0, representative = "mean_abs")

  expect_s4_class(msk, "LogicalNeuroVol")
  expect_equal(dim(msk), dims[1:3])
  expect_true(all(as.array(msk)[5:9, 5:9, 5:9]))
})

test_that("mapped and file-backed clip_level paths match dense results", {
  set.seed(1)
  fixture <- make_mask_backend_fixture()
  on.exit(unlink(fixture$path), add = TRUE)

  dense_scalar <- clip_level(fixture$dense, representative = "median", gradual = FALSE)
  dense_grad <- as.array(clip_level(fixture$dense, representative = "median", gradual = TRUE))

  expect_equal(
    clip_level(fixture$mmap, representative = "median", gradual = FALSE),
    dense_scalar
  )
  expect_equal(
    clip_level(fixture$filebacked, representative = "median", gradual = FALSE),
    dense_scalar
  )
  expect_equal(
    as.array(clip_level(fixture$mmap, representative = "median", gradual = TRUE)),
    dense_grad
  )
  expect_equal(
    as.array(clip_level(fixture$filebacked, representative = "median", gradual = TRUE)),
    dense_grad
  )
})

test_that("mapped and file-backed automask/apply_mask paths match dense results", {
  set.seed(2)
  fixture <- make_mask_backend_fixture()
  on.exit(unlink(fixture$path), add = TRUE)

  dense_mask <- automask(fixture$dense, representative = "mean_abs", gradual = FALSE, peels = 0)
  dense_mask_arr <- as.array(dense_mask)
  dense_applied <- as.array(apply_mask(fixture$dense, dense_mask))

  mmap_mask <- automask(fixture$mmap, representative = "mean_abs", gradual = FALSE, peels = 0)
  fb_mask <- automask(fixture$filebacked, representative = "mean_abs", gradual = FALSE, peels = 0)

  expect_equal(as.array(mmap_mask), dense_mask_arr)
  expect_equal(as.array(fb_mask), dense_mask_arr)
  expect_equal(
    unname(as.vector(as.array(apply_mask(fixture$mmap, dense_mask)))),
    unname(as.vector(dense_applied))
  )
  expect_equal(
    unname(as.vector(as.array(apply_mask(fixture$filebacked, dense_mask)))),
    unname(as.vector(dense_applied))
  )
})

test_that("apply_mask composes with automask", {
  dims <- c(11, 11, 11)
  sp <- NeuroSpace(dims)
  arr <- array(0, dims)
  arr[4:10, 4:10, 4:10] <- array(seq(50, 90, length.out = 7^3), c(7, 7, 7))
  vol <- NeuroVol(arr, sp)

  msk <- automask(vol, gradual = FALSE, peels = 0)
  masked <- apply_mask(vol, msk)

  expect_s4_class(masked, "DenseNeuroVol")
  expect_true(all(as.array(masked)[!as.array(msk)] == 0))
  expect_true(any(as.array(masked)[as.array(msk)] > 0))
})
