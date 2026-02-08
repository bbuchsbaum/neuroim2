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
