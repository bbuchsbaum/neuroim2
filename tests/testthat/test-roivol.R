
test_that("can create a spherical roi", {
  sp1 <- NeuroSpace(c(10,10,10), c(1,1,1))
  sphere <- spherical_roi(sp1, c(5,5,5), 3)
  expect_that(sphere, is_a("ROIVol"))
})

test_that("can create a cuboid roi", {
  sp1 <- NeuroSpace(c(10,10,10), c(1,1,1))
  cube <- cuboid_roi(sp1, c(5,5,5), 3)
  expect_that(cube, is_a("ROIVol"))
})

test_that("can convert roi coordinates to indices", {
  sp1 <- NeuroSpace(c(10,10,10), c(1,1,1))
  cube <- cuboid_roi(sp1, c(5,5,5), 3)
  idx <- grid_to_index(space(cube), coords(cube))
  vox <- index_to_grid(space(cube), idx)
  expect_equivalent(vox, coords(cube))
})

test_that("can add two ROIVols", {
  sp1 <- NeuroSpace(c(10,10,10), c(1,1,1))
  cube <- cuboid_roi(sp1, c(5,5,5), 3, fill=3)
  cube2 <- cuboid_roi(sp1, c(5,5,5), 3, fill=6)
  cube3 <- cube + cube2
  expect_true(all(cube3@.Data == 9))
})

