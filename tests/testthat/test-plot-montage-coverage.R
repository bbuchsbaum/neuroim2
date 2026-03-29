context("plot_montage coverage")

library(neuroim2)

test_that("plot_montage produces a ggplot for a DenseNeuroVol", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(1000), c(10, 10, 10)), sp)
  p <- plot_montage(vol)
  expect_s3_class(p, "gg")
})

test_that("plot_montage works with explicit zlevels", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(1000), c(10, 10, 10)), sp)
  p <- plot_montage(vol, zlevels = c(2L, 5L, 8L))
  expect_s3_class(p, "gg")
})

test_that("plot_montage works with range='data'", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(1000), c(10, 10, 10)), sp)
  p <- plot_montage(vol, range = "data")
  expect_s3_class(p, "gg")
})

test_that("plot_montage works with ncol and title/subtitle/caption", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(1000), c(10, 10, 10)), sp)
  p <- plot_montage(vol, ncol = 3L, title = "T", subtitle = "S", caption = "C")
  expect_s3_class(p, "gg")
})

test_that("plot_montage works along different axes", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(1000), c(10, 10, 10)), sp)
  p1 <- plot_montage(vol, along = 1L, zlevels = c(3L, 6L))
  expect_s3_class(p1, "gg")
  p2 <- plot_montage(vol, along = 2L, zlevels = c(3L, 6L))
  expect_s3_class(p2, "gg")
})

test_that("plot_montage works with downsample > 1", {
  sp <- NeuroSpace(c(20L, 20L, 20L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(8000), c(20, 20, 20)), sp)
  p <- plot_montage(vol, downsample = 2L, zlevels = c(5L, 10L))
  expect_s3_class(p, "gg")
})
