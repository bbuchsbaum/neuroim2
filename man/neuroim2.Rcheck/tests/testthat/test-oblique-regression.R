# Regression tests for oblique/sheared affine handling.
# Motivated by bugs fixed in v0.9.1 (plot memory blowup) and v0.10.0
# (downsample affine drift).

test_that("downsample preserves world-coordinate center for oblique affine", {
  skip_on_cran()
  aff <- matrix(c(2, 0.3, 0, -80,
                   0.3, 2, 0, -120,
                   0, 0, 2.5, -60,
                   0, 0, 0, 1), 4, 4, byrow = TRUE)
  sp <- NeuroSpace(c(32, 32, 16), trans = aff)
  vol <- DenseNeuroVol(rnorm(32 * 32 * 16), sp)

  vol2 <- downsample(vol, factor = 0.5)
  # Center voxel in world coords should be approximately preserved
  c1 <- grid_to_coord(space(vol), matrix(dim(vol) / 2, nrow = 1))
  c2 <- grid_to_coord(space(vol2), matrix(dim(vol2) / 2, nrow = 1))
  expect_equal(as.numeric(c1), as.numeric(c2), tolerance = max(spacing(space(vol))) * 2)
})

test_that("downsample of oblique vol has correct output dimensions", {
  skip_on_cran()
  aff <- matrix(c(1.5, 0.2, 0, -50,
                   0.0, 2.0, 0.1, -80,
                   0, 0, 2.5, -40,
                   0, 0, 0, 1), 4, 4, byrow = TRUE)
  sp <- NeuroSpace(c(20, 20, 10), trans = aff)
  vol <- DenseNeuroVol(rnorm(20 * 20 * 10), sp)

  vol2 <- downsample(vol, factor = 0.5)
  expect_equal(dim(vol2), as.integer(c(10, 10, 5)))
})

test_that("resample preserves data for identity transform", {
  skip_on_cran()
  vol <- make_vol(c(10, 10, 10), spacing = c(2, 2, 2))
  target_sp <- space(vol)

  vol2 <- resample(vol, target_sp)
  # Compare numeric values (ignoring object attributes from resample)
  expect_equal(as.numeric(vol@.Data), as.numeric(vol2@.Data), tolerance = 1e-4)
})

test_that("resample with oblique target space produces valid output", {
  skip_on_cran()
  vol <- make_vol(c(16, 16, 16), spacing = c(2, 2, 2))
  # Oblique target
  aff <- matrix(c(2, 0.1, 0, -16,
                   0.1, 2, 0, -16,
                   0, 0, 2, -16,
                   0, 0, 0, 1), 4, 4, byrow = TRUE)
  target_sp <- NeuroSpace(c(16, 16, 16), trans = aff)

  vol2 <- resample(vol, target_sp)
  expect_equal(dim(vol2), c(16L, 16L, 16L))
  expect_false(all(vol2@.Data == 0))
})

test_that("deoblique produces axis-aligned output", {
  skip_on_cran()
  aff <- matrix(c(2, 0.5, 0, -80,
                   0.0, 2, 0, -120,
                   0, 0, 2.5, -60,
                   0, 0, 0, 1), 4, 4, byrow = TRUE)
  sp <- NeuroSpace(c(32, 32, 16), trans = aff)
  vol <- DenseNeuroVol(rnorm(32 * 32 * 16), sp)

  dvol <- deoblique(vol)
  # After deoblique, off-diagonal elements should be ~0
  tx <- trans(space(dvol))
  off_diag <- c(tx[1,2], tx[1,3], tx[2,1], tx[2,3], tx[3,1], tx[3,2])
  expect_true(all(abs(off_diag) < 0.01))
})

test_that("deoblique(NeuroSpace) returns axis-aligned space", {
  aff <- matrix(c(2, 0.3, 0, -80,
                   0.3, 2, 0, -120,
                   0, 0, 2.5, -60,
                   0, 0, 0, 1), 4, 4, byrow = TRUE)
  sp <- NeuroSpace(c(32, 32, 16), trans = aff)

  dsp <- deoblique(sp)
  tx <- trans(dsp)
  off_diag <- c(tx[1,2], tx[1,3], tx[2,1], tx[2,3], tx[3,1], tx[3,2])
  expect_true(all(abs(off_diag) < 0.01))
})
