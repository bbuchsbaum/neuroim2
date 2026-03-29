context("plot_ortho coverage")

library(neuroim2)

test_that("plot_ortho returns list of ggplots invisibly", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(1000), c(10, 10, 10)), sp)
  # plot_ortho prints to grid and returns invisibly; capture output
  pdf(file = NULL)
  result <- plot_ortho(vol)
  dev.off()
  expect_true(is.list(result))
  expect_true(all(c("axial", "coronal", "sagittal") %in% names(result)))
})

test_that("plot_ortho panels are ggplot objects", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(1000), c(10, 10, 10)), sp)
  pdf(file = NULL)
  result <- plot_ortho(vol)
  dev.off()
  expect_s3_class(result$axial,    "gg")
  expect_s3_class(result$coronal,  "gg")
  expect_s3_class(result$sagittal, "gg")
})

test_that("plot_ortho works with explicit coord", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(1000), c(10, 10, 10)), sp)
  pdf(file = NULL)
  result <- plot_ortho(vol, coord = c(5L, 5L, 5L))
  dev.off()
  expect_true(is.list(result))
})

test_that("plot_ortho works with crosshair=FALSE and annotate=FALSE", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(1000), c(10, 10, 10)), sp)
  pdf(file = NULL)
  result <- plot_ortho(vol, crosshair = FALSE, annotate = FALSE)
  dev.off()
  expect_s3_class(result$axial, "gg")
})

test_that("plot_ortho works with range='data'", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(1000), c(10, 10, 10)), sp)
  pdf(file = NULL)
  result <- plot_ortho(vol, range = "data")
  dev.off()
  expect_s3_class(result$axial, "gg")
})

test_that("plot_ortho works with downsample > 1", {
  sp <- NeuroSpace(c(20L, 20L, 20L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(8000), c(20, 20, 20)), sp)
  pdf(file = NULL)
  result <- plot_ortho(vol, downsample = 2L)
  dev.off()
  expect_true(is.list(result))
})

test_that("plot_ortho unit='mm' path runs without error", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(2, 2, 2))
  vol <- DenseNeuroVol(array(rnorm(1000), c(10, 10, 10)), sp)
  pdf(file = NULL)
  # coord in mm — coord_to_grid is available in neuroim2
  result <- plot_ortho(vol, coord = c(5, 5, 5), unit = "mm")
  dev.off()
  expect_true(is.list(result))
})
