library(testthat)
library(neuroim2)

# ---- gaussian_blur ----

test_that("gaussian_blur without mask uses full volume", {
  sp  <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(1000), c(10, 10, 10)), sp)
  blurred <- gaussian_blur(vol, sigma = 2, window = 1)
  expect_s4_class(blurred, "NeuroVol")
  expect_equal(dim(blurred), dim(vol))
})

test_that("gaussian_blur smooths variance", {
  sp  <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  set.seed(42)
  vol <- DenseNeuroVol(array(rnorm(1000), c(10, 10, 10)), sp)
  blurred <- gaussian_blur(vol, sigma = 2, window = 1)
  expect_true(var(as.numeric(blurred@.Data)) < var(as.numeric(vol@.Data)))
})

test_that("gaussian_blur with LogicalNeuroVol mask", {
  sp   <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  vol  <- DenseNeuroVol(array(rnorm(1000), c(10, 10, 10)), sp)
  mask <- LogicalNeuroVol(array(TRUE, c(10, 10, 10)), sp)
  blurred <- gaussian_blur(vol, mask = mask, sigma = 2, window = 1)
  expect_s4_class(blurred, "NeuroVol")
  expect_equal(dim(blurred), dim(vol))
})

test_that("gaussian_blur with larger window", {
  sp  <- NeuroSpace(c(12L, 12L, 12L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(12^3), c(12, 12, 12)), sp)
  blurred <- gaussian_blur(vol, sigma = 1.5, window = 2)
  expect_s4_class(blurred, "NeuroVol")
})

test_that("gaussian_blur rejects invalid vol", {
  expect_error(gaussian_blur(matrix(1:9, 3, 3)), class = "error")
})

test_that("gaussian_blur rejects window < 1", {
  sp  <- NeuroSpace(c(8L, 8L, 8L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(512), c(8, 8, 8)), sp)
  expect_error(gaussian_blur(vol, sigma = 2, window = 0), class = "error")
})

test_that("gaussian_blur rejects sigma <= 0", {
  sp  <- NeuroSpace(c(8L, 8L, 8L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(512), c(8, 8, 8)), sp)
  expect_error(gaussian_blur(vol, sigma = -1, window = 1), class = "error")
})

test_that("gaussian_blur rejects invalid mask type", {
  sp   <- NeuroSpace(c(8L, 8L, 8L), c(1, 1, 1))
  vol  <- DenseNeuroVol(array(rnorm(512), c(8, 8, 8)), sp)
  expect_error(gaussian_blur(vol, mask = matrix(1:9, 3, 3)), class = "error")
})

# ---- guided_filter ----

test_that("guided_filter returns DenseNeuroVol with correct dims", {
  sp  <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(1000), c(10, 10, 10)), sp)
  out <- guided_filter(vol, radius = 2, epsilon = 0.1)
  expect_s4_class(out, "NeuroVol")
  expect_equal(dim(out), dim(vol))
})

test_that("guided_filter with default parameters", {
  sp  <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(1000), c(10, 10, 10)), sp)
  out <- guided_filter(vol)
  expect_s4_class(out, "NeuroVol")
})

test_that("guided_filter rejects invalid vol", {
  expect_error(guided_filter(array(1:27, c(3, 3, 3))), class = "error")
})

test_that("guided_filter rejects radius < 1", {
  sp  <- NeuroSpace(c(8L, 8L, 8L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(512), c(8, 8, 8)), sp)
  expect_error(guided_filter(vol, radius = 0), class = "error")
})

test_that("guided_filter rejects epsilon <= 0", {
  sp  <- NeuroSpace(c(8L, 8L, 8L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(512), c(8, 8, 8)), sp)
  expect_error(guided_filter(vol, epsilon = 0), class = "error")
})

# ---- bilateral_filter ----

test_that("bilateral_filter with explicit mask", {
  sp   <- NeuroSpace(c(8L, 8L, 8L), c(1, 1, 1))
  vol  <- DenseNeuroVol(array(rnorm(512), c(8, 8, 8)), sp)
  mask <- LogicalNeuroVol(array(TRUE, c(8, 8, 8)), sp)
  out  <- bilateral_filter(vol, mask = mask, spatial_sigma = 2, intensity_sigma = 1)
  expect_s4_class(out, "NeuroVol")
  expect_equal(dim(out), dim(vol))
})

test_that("bilateral_filter rejects window < 1", {
  sp  <- NeuroSpace(c(8L, 8L, 8L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(512), c(8, 8, 8)), sp)
  expect_error(bilateral_filter(vol, window = 0), class = "error")
})

test_that("bilateral_filter rejects spatial_sigma <= 0", {
  sp  <- NeuroSpace(c(8L, 8L, 8L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(512), c(8, 8, 8)), sp)
  expect_error(bilateral_filter(vol, spatial_sigma = 0), class = "error")
})

test_that("bilateral_filter rejects intensity_sigma <= 0", {
  sp  <- NeuroSpace(c(8L, 8L, 8L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(512), c(8, 8, 8)), sp)
  expect_error(bilateral_filter(vol, intensity_sigma = -1), class = "error")
})

test_that("bilateral_filter rejects invalid mask type", {
  sp  <- NeuroSpace(c(8L, 8L, 8L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(512), c(8, 8, 8)), sp)
  expect_error(bilateral_filter(vol, mask = array(TRUE, c(8, 8, 8))), class = "error")
})

# ---- bilateral_filter_vec (internal function via :::) ----

test_that("bilateral_filter_vec without mask returns NeuroVec", {
  sp  <- NeuroSpace(c(8L, 8L, 8L, 4L), c(1, 1, 1))
  vec <- DenseNeuroVec(array(rnorm(8 * 8 * 8 * 4), c(8, 8, 8, 4)), sp)
  out <- neuroim2:::bilateral_filter_vec(vec, spatial_sigma = 1, intensity_sigma = 1, window = 1)
  expect_s4_class(out, "NeuroVec")
  expect_equal(dim(out)[1:3], c(8L, 8L, 8L))
  expect_equal(dim(out)[4], 4L)
})

test_that("bilateral_filter_vec with mask", {
  sp   <- NeuroSpace(c(8L, 8L, 8L, 4L), c(1, 1, 1))
  vec  <- DenseNeuroVec(array(rnorm(8 * 8 * 8 * 4), c(8, 8, 8, 4)), sp)
  msp  <- NeuroSpace(c(8L, 8L, 8L), c(1, 1, 1))
  mask <- LogicalNeuroVol(array(TRUE, c(8, 8, 8)), msp)
  out  <- neuroim2:::bilateral_filter_vec(vec, mask = mask, spatial_sigma = 1, intensity_sigma = 1)
  expect_s4_class(out, "NeuroVec")
})

test_that("bilateral_filter_vec rejects non-NeuroVec", {
  expect_error(neuroim2:::bilateral_filter_vec(array(1:27, c(3, 3, 3))), class = "error")
})

test_that("bilateral_filter_vec rejects window < 1", {
  sp  <- NeuroSpace(c(6L, 6L, 6L, 3L), c(1, 1, 1))
  vec <- DenseNeuroVec(array(rnorm(6^3 * 3), c(6, 6, 6, 3)), sp)
  expect_error(neuroim2:::bilateral_filter_vec(vec, window = 0), class = "error")
})

test_that("bilateral_filter_vec rejects invalid mask type", {
  sp  <- NeuroSpace(c(6L, 6L, 6L, 3L), c(1, 1, 1))
  vec <- DenseNeuroVec(array(rnorm(6^3 * 3), c(6, 6, 6, 3)), sp)
  expect_error(neuroim2:::bilateral_filter_vec(vec, mask = array(TRUE, c(6, 6, 6))), class = "error")
})

# ---- bilateral_filter_4d ----

test_that("bilateral_filter_4d without mask returns DenseNeuroVec", {
  sp  <- NeuroSpace(c(6L, 6L, 6L, 4L), c(1, 1, 1))
  vec <- DenseNeuroVec(array(rnorm(6^3 * 4), c(6, 6, 6, 4)), sp)
  out <- bilateral_filter_4d(vec, spatial_sigma = 1, intensity_sigma = 1,
                              temporal_sigma = 1, spatial_window = 1, temporal_window = 1)
  expect_s4_class(out, "DenseNeuroVec")
  expect_equal(dim(out), dim(vec))
})

test_that("bilateral_filter_4d with mask", {
  sp   <- NeuroSpace(c(6L, 6L, 6L, 4L), c(1, 1, 1))
  vec  <- DenseNeuroVec(array(rnorm(6^3 * 4), c(6, 6, 6, 4)), sp)
  msp  <- NeuroSpace(c(6L, 6L, 6L), c(1, 1, 1))
  mask <- LogicalNeuroVol(array(TRUE, c(6, 6, 6)), msp)
  out  <- bilateral_filter_4d(vec, mask = mask, spatial_sigma = 1, intensity_sigma = 1,
                               temporal_sigma = 1, spatial_window = 1, temporal_window = 1)
  expect_s4_class(out, "DenseNeuroVec")
})

test_that("bilateral_filter_4d rejects non-NeuroVec", {
  expect_error(bilateral_filter_4d(array(1:24, c(2, 3, 4, 1))), class = "error")
})

test_that("bilateral_filter_4d rejects spatial_window < 1", {
  sp  <- NeuroSpace(c(6L, 6L, 6L, 3L), c(1, 1, 1))
  vec <- DenseNeuroVec(array(rnorm(6^3 * 3), c(6, 6, 6, 3)), sp)
  expect_error(bilateral_filter_4d(vec, spatial_window = 0), class = "error")
})

test_that("bilateral_filter_4d rejects temporal_window < 0", {
  sp  <- NeuroSpace(c(6L, 6L, 6L, 3L), c(1, 1, 1))
  vec <- DenseNeuroVec(array(rnorm(6^3 * 3), c(6, 6, 6, 3)), sp)
  expect_error(bilateral_filter_4d(vec, temporal_window = -1), class = "error")
})

test_that("bilateral_filter_4d rejects spatial_sigma <= 0", {
  sp  <- NeuroSpace(c(6L, 6L, 6L, 3L), c(1, 1, 1))
  vec <- DenseNeuroVec(array(rnorm(6^3 * 3), c(6, 6, 6, 3)), sp)
  expect_error(bilateral_filter_4d(vec, spatial_sigma = 0), class = "error")
})

test_that("bilateral_filter_4d rejects temporal_sigma <= 0", {
  sp  <- NeuroSpace(c(6L, 6L, 6L, 3L), c(1, 1, 1))
  vec <- DenseNeuroVec(array(rnorm(6^3 * 3), c(6, 6, 6, 3)), sp)
  expect_error(bilateral_filter_4d(vec, temporal_sigma = 0), class = "error")
})

test_that("bilateral_filter_4d rejects intensity_sigma <= 0", {
  sp  <- NeuroSpace(c(6L, 6L, 6L, 3L), c(1, 1, 1))
  vec <- DenseNeuroVec(array(rnorm(6^3 * 3), c(6, 6, 6, 3)), sp)
  expect_error(bilateral_filter_4d(vec, intensity_sigma = -0.1), class = "error")
})

test_that("bilateral_filter_4d rejects temporal_spacing <= 0", {
  sp  <- NeuroSpace(c(6L, 6L, 6L, 3L), c(1, 1, 1))
  vec <- DenseNeuroVec(array(rnorm(6^3 * 3), c(6, 6, 6, 3)), sp)
  expect_error(bilateral_filter_4d(vec, temporal_spacing = 0), class = "error")
})

# ---- laplace_enhance ----

test_that("laplace_enhance without mask returns NeuroVol", {
  sp  <- NeuroSpace(c(8L, 8L, 8L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(512), c(8, 8, 8)), sp)
  out <- laplace_enhance(vol, k = 1, patch_size = 3, search_radius = 1, h = 0.5)
  expect_s4_class(out, "NeuroVol")
  expect_equal(dim(out), dim(vol))
})

test_that("laplace_enhance with explicit LogicalNeuroVol mask", {
  sp   <- NeuroSpace(c(8L, 8L, 8L), c(1, 1, 1))
  vol  <- DenseNeuroVol(array(rnorm(512), c(8, 8, 8)), sp)
  mask <- LogicalNeuroVol(array(TRUE, c(8, 8, 8)), sp)
  out  <- laplace_enhance(vol, mask = mask, k = 1, patch_size = 3, search_radius = 1)
  expect_s4_class(out, "NeuroVol")
})

test_that("laplace_enhance rejects non-NeuroVol", {
  expect_error(laplace_enhance(array(1:27, c(3, 3, 3))), class = "error")
})

test_that("laplace_enhance rejects k < 1", {
  sp  <- NeuroSpace(c(8L, 8L, 8L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(512), c(8, 8, 8)), sp)
  expect_error(laplace_enhance(vol, k = 0), class = "error")
})

test_that("laplace_enhance rejects even patch_size", {
  sp  <- NeuroSpace(c(8L, 8L, 8L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(512), c(8, 8, 8)), sp)
  expect_error(laplace_enhance(vol, patch_size = 4), class = "error")
})

test_that("laplace_enhance rejects patch_size < 3", {
  sp  <- NeuroSpace(c(8L, 8L, 8L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(512), c(8, 8, 8)), sp)
  expect_error(laplace_enhance(vol, patch_size = 1), class = "error")
})

test_that("laplace_enhance rejects search_radius < 1", {
  sp  <- NeuroSpace(c(8L, 8L, 8L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(512), c(8, 8, 8)), sp)
  expect_error(laplace_enhance(vol, search_radius = 0), class = "error")
})

test_that("laplace_enhance rejects h <= 0", {
  sp  <- NeuroSpace(c(8L, 8L, 8L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(512), c(8, 8, 8)), sp)
  expect_error(laplace_enhance(vol, h = 0), class = "error")
})

test_that("laplace_enhance rejects non-LogicalNeuroVol mask", {
  sp   <- NeuroSpace(c(8L, 8L, 8L), c(1, 1, 1))
  vol  <- DenseNeuroVol(array(rnorm(512), c(8, 8, 8)), sp)
  expect_error(laplace_enhance(vol, mask = array(TRUE, c(8, 8, 8))), class = "error")
})
