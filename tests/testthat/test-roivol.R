

context("roivol")


test_that("can create a spherical roi", {
  sp1 <- NeuroSpace(c(10,10,10), c(1,1,1))
  sphere <- spherical_roi(sp1, c(5,5,5), 3)
  expect_that(sphere, is_a("ROIVol"))
})

test_that("spherical roi with too small radius causes error", {
  sp1 <- NeuroSpace(c(10,10,10), c(1,1,1))
  expect_error(spherical_roi(sp1, c(5,5,5), .5))

})

test_that("can create a cuboid roi", {
  sp1 <- NeuroSpace(c(10,10,10), c(1,1,1))
  cube <- cuboid_roi(sp1, c(5,5,5), 3)
  expect_that(cube, is_a("ROIVol"))
})

test_that("cuboid roi must have three elements", {
  sp1 <- NeuroSpace(c(10,10,10), c(1,1,1))
  cube <- expect_error(cuboid_roi(sp1, c(5,5,5,6), 3))
})

test_that("cuboid roi must have positive surround", {
  sp1 <- NeuroSpace(c(10,10,10), c(1,1,1))
  cube <- expect_error(cuboid_roi(sp1, c(5,5,5), -1))
})

test_that("square roi must have three elements", {
  sp1 <- NeuroSpace(c(10,10,10), c(1,1,1))
  cube <- expect_error(cuboid_roi(sp1, c(5,5,5,6), 3))
})

test_that("square roi must have positive surround", {
  sp1 <- NeuroSpace(c(10,10,10), c(1,1,1))
  cube <- expect_error(square_roi(sp1, c(5,5,5), -1))
})




test_that("can create a square roi", {
  sp1 <- NeuroSpace(c(10,10,10), c(1,1,1))
  square <- square_roi(sp1, c(5,5,5), 1)
  expect_that(square, is_a("ROIVol"))
  expect_equal(nrow(coords(square)), 9)
})

test_that("can convert roi coordinates to indices", {
  sp1 <- NeuroSpace(c(10,10,10), c(1,1,1))
  cube <- cuboid_roi(sp1, c(5,5,5), 3)
  idx <- grid_to_index(space(cube), coords(cube))
  vox <- index_to_grid(space(cube), idx)
  expect_equivalent(vox, coords(cube))
})
test_that("can convert ROIVol to DenseNeuroVol", {
  sp1 <- NeuroSpace(c(10,10,10), c(1,1,1))
  cube <- cuboid_roi(sp1, c(5,5,5), 3)
  vol <- as.dense(cube)
  expect_equal(sum(vol), sum(cube))
})

test_that("can add two ROIVols", {
  sp1 <- NeuroSpace(c(10,10,10), c(1,1,1))
  cube <- cuboid_roi(sp1, c(5,5,5), 3, fill=3)
  cube2 <- cuboid_roi(sp1, c(5,5,5), 3, fill=6)
  cube3 <- cube + cube2
  expect_true(all(cube3@.Data == 9))
})

test_that("can construct an ROIVec", {
  sp1 <- NeuroSpace(c(10,10,10,4), c(1,1,1))
  cube <- cuboid_roi(drop_dim(sp1), c(5,5,5), 3, fill=3)
  vec <- ROIVec(sp1, coords(cube), matrix(1, 5, nrow(coords(cube))))
  expect_true(!is.null(vec))
})

test_that("can convert ROIVec to matrix", {
  sp1 <- NeuroSpace(c(10,10,10,4), c(1,1,1))
  cube <- cuboid_roi(drop_dim(sp1), c(5,5,5), 3, fill=3)
  m <- matrix(1, 5, nrow(coords(cube)))
  vec <- ROIVec(sp1, coords(cube), matrix(1, 5, nrow(coords(cube))))
  expect_true(inherits(as(vec, "matrix"), "matrix"))
})


test_that("square_roi keeps center voxel when nonzero=TRUE even if center is zero", {
  sp <- NeuroSpace(c(3,3,3), c(1,1,1))
  arr <- array(0, dim=c(3,3,3))
  arr[1,1,1] <- 5  # one nonzero voxel
  vol <- NeuroVol(arr, sp)
  centroid <- c(2,2,2)

  roi <- square_roi(vol, centroid, surround=1, nonzero=TRUE)

  # center voxel is retained
  expect_true(any(apply(coords(roi), 1, function(r) all(r == centroid))))
  expect_length(roi@center_index, 1)
  expect_false(is.na(roi@center_index))
  expect_equal(roi@parent_index, grid_to_index(vol, matrix(centroid, ncol=3)))
  # zero-valued center should remain despite nonzero filter
  expect_true(any(values(roi) == 0))
})

test_that("cuboid_roi keeps center voxel when nonzero=TRUE even if center is zero", {
  sp <- NeuroSpace(c(3,3,3), c(1,1,1))
  arr <- array(0, dim=c(3,3,3))
  arr[1,1,1] <- 5  # one nonzero voxel
  vol <- NeuroVol(arr, sp)
  centroid <- c(2,2,2)

  roi <- cuboid_roi(vol, centroid, surround=1, nonzero=TRUE)

  expect_true(any(apply(coords(roi), 1, function(r) all(r == centroid))))
  expect_length(roi@center_index, 1)
  expect_false(is.na(roi@center_index))
  expect_equal(roi@parent_index, grid_to_index(vol, matrix(centroid, ncol=3)))
  expect_true(any(values(roi) == 0))
})


