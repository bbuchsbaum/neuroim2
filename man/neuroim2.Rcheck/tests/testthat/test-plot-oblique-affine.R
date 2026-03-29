library(testthat)

test_that("plot(NeuroVol) is safe for oblique affine spaces", {
  dims <- c(24L, 28L, 12L)
  trans <- diag(4)
  trans[1, 2] <- 0.25
  trans[2, 1] <- -0.15
  sp <- neuroim2::NeuroSpace(dims, trans = trans)
  vol <- neuroim2::NeuroVol(array(stats::rnorm(prod(dims)), dim = dims), sp)

  p <- neuroim2::plot(vol, zlevels = c(3L, 6L, 9L))
  expect_s3_class(p, "ggplot")

  tf <- tempfile(fileext = ".png")
  grDevices::png(tf)
  on.exit({
    grDevices::dev.off()
    unlink(tf)
  }, add = TRUE)

  expect_warning(print(p), NA)
})

test_that("plot(NeuroSlice) is safe for oblique affine spaces", {
  dims <- c(24L, 28L, 12L)
  trans <- diag(4)
  trans[1, 2] <- 0.2
  trans[2, 1] <- -0.1
  sp <- neuroim2::NeuroSpace(dims, trans = trans)
  vol <- neuroim2::NeuroVol(array(stats::rnorm(prod(dims)), dim = dims), sp)
  slc <- neuroim2::slice(vol, zlevel = 6L, along = 3L)

  p <- neuroim2::plot(slc)
  expect_s3_class(p, "ggplot")

  tf <- tempfile(fileext = ".png")
  grDevices::png(tf)
  on.exit({
    grDevices::dev.off()
    unlink(tf)
  }, add = TRUE)

  expect_warning(print(p), NA)
})
