library(testthat)

test_that("deoblique() builds AFNI-style target space from oblique affine", {
  tx <- diag(c(2, 3, 4, 1))
  tx[1, 2] <- 0.25
  tx[2, 1] <- -0.1
  sp <- neuroim2::NeuroSpace(c(18, 12, 10), spacing = c(2, 3, 4), trans = tx)

  expect_gt(max(neuroim2::obliquity(neuroim2::trans(sp))), 0)

  sp_deob <- neuroim2::deoblique(sp)
  expect_s4_class(sp_deob, "NeuroSpace")
  expect_equal(diag(neuroim2::trans(sp_deob))[1:3], rep(2, 3))
  expect_equal(max(neuroim2::obliquity(neuroim2::trans(sp_deob))), 0, tolerance = 1e-10)
})

test_that("deoblique() supports newgrid and gridset controls", {
  tx <- diag(c(2, 3, 4, 1))
  tx[1, 3] <- 0.2
  sp <- neuroim2::NeuroSpace(c(16, 14, 8), spacing = c(2, 3, 4), trans = tx)

  sp_newgrid <- neuroim2::deoblique(sp, newgrid = 1.5)
  expect_equal(diag(neuroim2::trans(sp_newgrid))[1:3], rep(1.5, 3))

  grid_sp <- neuroim2::NeuroSpace(c(20, 20, 20), spacing = c(1, 1, 1), origin = c(-5, -6, -7))
  sp_gridset <- neuroim2::deoblique(sp, gridset = grid_sp)
  expect_equal(dim(sp_gridset), dim(grid_sp))
  expect_equal(neuroim2::trans(sp_gridset), neuroim2::trans(grid_sp))

  grid_sp_2d <- neuroim2::NeuroSpace(c(20, 20), spacing = c(1, 1))
  expect_error(neuroim2::deoblique(sp, gridset = grid_sp_2d), "3D")

  expect_error(neuroim2::deoblique(sp, gridset = grid_sp, newgrid = 2), "mutually exclusive")
})

test_that("deoblique() resamples NeuroVol to deobliqued target", {
  tx <- diag(c(2, 3, 4, 1))
  tx[1, 2] <- 0.3
  sp <- neuroim2::NeuroSpace(c(14, 10, 8), spacing = c(2, 3, 4), trans = tx)
  vol <- neuroim2::NeuroVol(array(stats::rnorm(prod(dim(sp))), dim = dim(sp)), sp)

  vol_deob <- neuroim2::deoblique(vol, newgrid = 2, method = "linear")
  expect_s4_class(vol_deob, "NeuroVol")
  expect_equal(diag(neuroim2::trans(neuroim2::space(vol_deob)))[1:3], rep(2, 3))
  expect_equal(max(neuroim2::obliquity(neuroim2::trans(neuroim2::space(vol_deob)))), 0, tolerance = 1e-10)
})

test_that("deoblique() errors on non-3D spaces", {
  sp2 <- neuroim2::NeuroSpace(c(20, 20), spacing = c(1, 1))
  expect_error(neuroim2::deoblique(sp2), "3D")
})
