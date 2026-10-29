library(testthat)
library(neuroim2)



# Test the resample(NeuroVol, NeuroVol) function
test_that("resample(NeuroVol, NeuroVol) works correctly", {
  skip_on_cran()
  source <- NeuroVol(array(rnorm(64 * 64 * 64), c(64, 64, 64)), NeuroSpace(c(64, 64, 64)))
  target <- NeuroVol(array(rnorm(64 * 64 * 64), c(64, 64, 64)), NeuroSpace(c(64, 64, 64)))
  interpolation <- 3L

  resampled_vol <- resample(source, target, interpolation)

  expect_s4_class(resampled_vol, "NeuroVol")
  expect_equal(dim(resampled_vol), dim(target))
  expect_equal(space(resampled_vol), space(target))
})

# Test the resample(NeuroVol, NeuroSpace) function
test_that("resample(NeuroVol, NeuroSpace) works correctly", {
  skip_on_cran()
  source <- NeuroVol(array(rnorm(64 * 64 * 64), c(64, 64, 64)), NeuroSpace(c(64, 64, 64)))
  target_space <- NeuroSpace(c(64, 64, 64))
  interpolation <- 3L

  resampled_vol <- resample(source, target_space, interpolation)

  expect_s4_class(resampled_vol, "NeuroVol")
  expect_equal(dim(resampled_vol), dim(target_space))
  expect_equal(space(resampled_vol), target_space)
})

test_that("resample(ClusteredNeuroVol, NeuroSpace) preserves labels and label map", {
  sp <- NeuroSpace(c(4, 4, 4))
  mask_data <- array(FALSE, dim = c(4, 4, 4))
  mask_data[1:2, 1:2, 1:2] <- TRUE
  mask_data[3:4, 3:4, 3:4] <- TRUE

  mask <- LogicalNeuroVol(mask_data, sp)
  clusters <- rep(c(1L, 2L), each = sum(mask_data) / 2)
  label_map <- list(regionA = 1L, regionB = 2L)
  cvol <- ClusteredNeuroVol(mask, clusters, label_map = label_map)

  target_space <- NeuroSpace(c(4, 4, 4), spacing = c(0.5, 0.5, 0.5))

  resampled <- resample(cvol, target_space)

  expect_s4_class(resampled, "ClusteredNeuroVol")
  expect_equal(space(resampled), target_space)
  expect_setequal(unique(resampled@clusters), c(1L, 2L))
  expect_equal(names(resampled@label_map), names(label_map))
})

test_that("resample(ClusteredNeuroVol, NeuroVol) forces nearest neighbor interpolation", {
  sp <- NeuroSpace(c(4, 4, 4))
  mask <- LogicalNeuroVol(array(TRUE, dim = c(4, 4, 4)), sp)
  clusters <- rep(1:2, each = 32)
  cvol <- ClusteredNeuroVol(mask, as.integer(clusters))

  target_vol <- NeuroVol(array(0, dim = c(4, 4, 4)), sp)

  expect_warning(
    res <- resample(cvol, target_vol, interpolation = 3L),
    "nearest-neighbor interpolation"
  )

  expect_s4_class(res, "ClusteredNeuroVol")
  expect_equal(space(res), space(target_vol))
  expect_setequal(unique(res@clusters), c(1L, 2L))
})

test_that("a volume resamples onto the spatial grid of a 4-D target", {
  skip_on_cran()
  source <- NeuroVol(array(rnorm(8 * 8 * 8), c(8, 8, 8)), NeuroSpace(c(8, 8, 8)))
  grid3 <- NeuroSpace(c(6, 6, 6), spacing = c(1.2, 1.2, 1.2), origin = c(1, 2, 3))
  run <- NeuroVec(array(0, c(6, 6, 6, 4)), add_dim(grid3, 4))
  expected <- resample(source, grid3, 1L)

  for (target in list(run, space(run))) {
    out <- resample(source, target, 1L)
    expect_s4_class(out, "NeuroVol")
    expect_equal(dim(out), c(6L, 6L, 6L))
    expect_equal(trans(space(out)), trans(grid3))
    expect_equal(as.array(out), as.array(expected))
  }
})

test_that("resample(ClusteredNeuroVol, <4-D target>) keeps labels on the spatial grid", {
  sp <- NeuroSpace(c(4, 4, 4))
  mask <- LogicalNeuroVol(array(TRUE, dim = c(4, 4, 4)), sp)
  cvol <- ClusteredNeuroVol(mask, rep(1:2, each = 32),
                            label_map = list(a = 1L, b = 2L))
  grid3 <- NeuroSpace(c(4, 4, 4), spacing = c(0.5, 0.5, 0.5))
  run <- NeuroVec(array(0, c(4, 4, 4, 3)), add_dim(grid3, 3))
  expected <- resample(cvol, grid3)

  for (target in list(run, space(run))) {
    out <- resample(cvol, target)
    expect_s4_class(out, "ClusteredNeuroVol")
    expect_equal(dim(out), c(4L, 4L, 4L))
    expect_equal(out@clusters, expected@clusters)
    expect_equal(names(out@label_map), c("a", "b"))
  }
})
