

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


# nonzero = TRUE means what it says: a zero-valued centre is dropped like any
# other zero voxel. It used to be forced back in so that parent_index could be
# read off the coordinate matrix; parent_index is now derived from the requested
# centroid directly, so the workaround is no longer needed -- and forcing a zero
# voxel into a "nonzero" ROI contradicted the argument's documentation.
test_that("square_roi drops a zero-valued center when nonzero=TRUE", {
  sp <- NeuroSpace(c(3,3,3), c(1,1,1))
  arr <- array(0, dim=c(3,3,3))
  arr[1,1,1] <- 5  # one nonzero voxel
  vol <- NeuroVol(arr, sp)
  centroid <- c(2,2,2)

  roi <- square_roi(vol, centroid, surround=1, nonzero=TRUE)

  expect_false(any(apply(coords(roi), 1, function(r) all(r == centroid))))
  expect_length(roi@center_index, 1)
  expect_true(is.na(roi@center_index))
  # parent_index still points at the requested centroid, filtering or not
  expect_equal(roi@parent_index, grid_to_index(vol, matrix(centroid, ncol=3)))
  # and no zero-valued voxel survives the filter
  expect_false(any(values(roi) == 0))
})

test_that("cuboid_roi drops a zero-valued center when nonzero=TRUE", {
  sp <- NeuroSpace(c(3,3,3), c(1,1,1))
  arr <- array(0, dim=c(3,3,3))
  arr[1,1,1] <- 5  # one nonzero voxel
  vol <- NeuroVol(arr, sp)
  centroid <- c(2,2,2)

  roi <- cuboid_roi(vol, centroid, surround=1, nonzero=TRUE)

  expect_false(any(apply(coords(roi), 1, function(r) all(r == centroid))))
  expect_length(roi@center_index, 1)
  expect_true(is.na(roi@center_index))
  expect_equal(roi@parent_index, grid_to_index(vol, matrix(centroid, ncol=3)))
  expect_false(any(values(roi) == 0))
})

test_that("the windowed ROI builders agree on centre and parent semantics", {
  sp <- NeuroSpace(c(10,10,10), c(1,1,1))
  arr <- array(1, dim=c(10,10,10))
  arr[5,5,5] <- 0
  vol <- NeuroVol(arr, sp)
  ctr <- c(5,5,5)
  pidx <- grid_to_index(vol, matrix(ctr, ncol=3))

  builders <- list(
    spherical = function(nz) spherical_roi(vol, ctr, 2.5, nonzero = nz),
    cuboid    = function(nz) cuboid_roi(vol, ctr, 2, nonzero = nz),
    square    = function(nz) square_roi(vol, ctr, 2, nonzero = nz)
  )

  for (nm in names(builders)) {
    keep <- builders[[nm]](FALSE)
    drop <- builders[[nm]](TRUE)

    # centre present when unfiltered, absent once its zero value is filtered
    expect_false(is.na(keep@center_index), info = nm)
    expect_true(is.na(drop@center_index), info = nm)
    expect_equal(nrow(drop@coords), nrow(keep@coords) - 1L, info = nm)

    # parent_index is the centroid either way
    expect_equal(keep@parent_index, pidx, info = nm)
    expect_equal(drop@parent_index, pidx, info = nm)

    # coords are integer, undecorated, and lexicographically ordered
    for (r in list(keep, drop)) {
      expect_identical(storage.mode(r@coords), "integer", info = nm)
      expect_null(dimnames(r@coords), info = nm)
      expect_identical(order(r@coords[,1], r@coords[,2], r@coords[,3]),
                       seq_len(nrow(r@coords)), info = nm)
    }
  }

  # and they all validate the centroid the same way
  for (nm in names(builders)) {
    f <- switch(nm, spherical = function(c) spherical_roi(vol, c, 2.5),
                    cuboid    = function(c) cuboid_roi(vol, c, 2),
                    square    = function(c) square_roi(vol, c, 2))
    expect_error(f(c(99, 5, 5)), "exceeds volume dimensions", info = nm)
    expect_error(f(c(0, 5, 5)), ">= 1", info = nm)
    expect_error(f(c(5, 5)), "length 3", info = nm)
  }
})


