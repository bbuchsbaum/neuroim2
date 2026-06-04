library(testthat)
library(neuroim2)

# Test the gaussian_blur function
test_that("gaussian_blur works correctly", {
  skip_on_cran()
  vol <- NeuroVol(array(rnorm(64 * 64 * 64), c(64, 64, 64)), NeuroSpace(c(64, 64, 64)))
  mask <- as.logical(NeuroVol(as.logical(round(runif(prod(dim(vol))))), space(vol)))
  sigma <- 2
  window <- 1

  blurred_vol <- gaussian_blur(vol, mask, sigma, window)

  expect_s4_class(blurred_vol, "NeuroVol")
  expect_equal(dim(blurred_vol), dim(vol))
  expect_true(all(blurred_vol[mask == 0] == 0))
})

# Masked blur must be insulated to the mask (GitHub issue #22)
test_that("gaussian_blur insulates the mask against out-of-mask NaN", {
  skip_on_cran()
  sp <- NeuroSpace(c(20, 20, 20), spacing = c(2, 2, 2))
  mk  <- array(FALSE, c(20, 20, 20)); mk[6:15, , ] <- TRUE
  arr <- array(NaN,   c(20, 20, 20)); arr[mk] <- 1            # NaN outside the mask

  out <- gaussian_blur(NeuroVol(arr, sp), LogicalNeuroVol(mk, sp),
                       sigma = 1.27, window = 4)
  out <- as.numeric(out); dim(out) <- c(20, 20, 20)

  # No in-mask voxel may be erased by out-of-mask NaN.
  expect_equal(sum(is.nan(out) & mk), 0L)
  # A constant in-mask map stays constant, edge included (renormalized).
  expect_equal(out[6, 10, 10], 1)
  expect_equal(out[10, 10, 10], 1)
})

test_that("gaussian_blur renormalizes mask edges (no zero-bias)", {
  skip_on_cran()
  sp <- NeuroSpace(c(20, 20, 20), spacing = c(2, 2, 2))
  mk  <- array(FALSE, c(20, 20, 20)); mk[6:15, , ] <- TRUE
  arr <- array(0, c(20, 20, 20)); arr[mk] <- 1               # finite 0 outside the mask

  out <- gaussian_blur(NeuroVol(arr, sp), LogicalNeuroVol(mk, sp),
                       sigma = 1.27, window = 4)
  out <- as.numeric(out); dim(out) <- c(20, 20, 20)

  expect_equal(out[6, 10, 10], 1)    # edge not pulled toward exterior 0
  expect_equal(out[10, 10, 10], 1)   # interior unchanged
})

test_that("gaussian_blur(normalize=TRUE) matches smooth-in-mask workaround", {
  skip_on_cran()
  sp <- NeuroSpace(c(20, 20, 20), spacing = c(2, 2, 2))
  mk  <- array(FALSE, c(20, 20, 20)); mk[6:15, , ] <- TRUE
  set.seed(1)
  dat <- array(NaN, c(20, 20, 20)); dat[mk] <- rnorm(sum(mk))
  full <- LogicalNeuroVol(array(TRUE, c(20, 20, 20)), sp)
  sigma <- 1.27; window <- 4L

  new_def <- as.numeric(gaussian_blur(NeuroVol(dat, sp), LogicalNeuroVol(mk, sp),
                                      sigma = sigma, window = window))
  m   <- is.finite(dat)
  num <- as.numeric(gaussian_blur(NeuroVol(ifelse(m, dat, 0), sp), full,
                                  sigma = sigma, window = window, normalize = FALSE))
  den <- as.numeric(gaussian_blur(NeuroVol(m * 1.0, sp), full,
                                  sigma = sigma, window = window, normalize = FALSE))
  sm  <- num / den; sm[!as.logical(m)] <- 0; sm[!is.finite(sm)] <- 0

  expect_equal(new_def[as.logical(mk)], sm[as.logical(mk)])
})

test_that("gaussian_blur(normalize=FALSE) preserves legacy full-kernel behavior", {
  skip_on_cran()
  sp <- NeuroSpace(c(20, 20, 20), spacing = c(2, 2, 2))
  mk  <- array(FALSE, c(20, 20, 20)); mk[6:15, , ] <- TRUE
  arr <- array(NaN, c(20, 20, 20)); arr[mk] <- 1

  out <- gaussian_blur(NeuroVol(arr, sp), LogicalNeuroVol(mk, sp),
                       sigma = 1.27, window = 4, normalize = FALSE)
  out <- as.numeric(out); dim(out) <- c(20, 20, 20)

  # Legacy path leaks out-of-mask NaN into the boundary shell.
  expect_gt(sum(is.nan(out) & mk), 0L)
})

test_that("gaussian_blur rejects invalid normalize", {
  sp  <- NeuroSpace(c(8L, 8L, 8L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(512), c(8, 8, 8)), sp)
  expect_error(gaussian_blur(vol, sigma = 2, window = 1, normalize = NA),
               class = "error")
  expect_error(gaussian_blur(vol, sigma = 2, window = 1, normalize = "yes"),
               class = "error")
})

# Test the guided_filter function
test_that("guided_filter works correctly", {
  skip_on_cran()
  vol <- NeuroVol(array(rnorm(32 * 32 * 32), c(32, 32, 32)), NeuroSpace(c(32, 32, 32)))
  radius <- 4
  epsilon <- 0.7^2

  filtered_vol <- guided_filter(vol, radius, epsilon)

  expect_s4_class(filtered_vol, "NeuroVol")
  expect_equal(dim(filtered_vol), dim(vol))
})

test_that("bilateral_filter handles missing mask", {
  skip_on_cran()
  vol <- NeuroVol(array(rnorm(20^3), c(20, 20, 20)), NeuroSpace(c(20, 20, 20)))

  filtered_vol <- bilateral_filter(vol, spatial_sigma = 2, intensity_sigma = 1, window = 1)

  expect_s4_class(filtered_vol, "NeuroVol")
  expect_equal(dim(filtered_vol), dim(vol))
  expect_identical(space(filtered_vol), space(vol))
})
