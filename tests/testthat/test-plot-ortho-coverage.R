context("plot_ortho coverage")

library(neuroim2)

test_that("plot_ortho returns one assembled figure (visible unless drawn)", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(1000), c(10, 10, 10)), sp)
  pdf(file = NULL)
  on.exit(dev.off(), add = TRUE)
  # ggplot idiom: without drawing, the figure is returned visibly so it
  # auto-prints; draw = TRUE prints it and returns it invisibly.
  expect_visible(plot_ortho(vol, draw = FALSE))
  expect_invisible(plot_ortho(vol, draw = TRUE))
  result <- plot_ortho(vol, draw = FALSE)
  expect_s3_class(result, "patchwork")
})

test_that("plot_ortho(assemble = FALSE) returns the named list of panels", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(1000), c(10, 10, 10)), sp)
  pdf(file = NULL)
  on.exit(dev.off(), add = TRUE)
  result <- plot_ortho(vol, assemble = FALSE)
  expect_true(is.list(result))
  expect_true(all(c("axial", "coronal", "sagittal") %in% names(result)))
})

test_that("plot_ortho panels are ggplot objects", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(1000), c(10, 10, 10)), sp)
  pdf(file = NULL)
  on.exit(dev.off(), add = TRUE)
  result <- plot_ortho(vol, assemble = FALSE)
  expect_s3_class(result$axial,    "gg")
  expect_s3_class(result$coronal,  "gg")
  expect_s3_class(result$sagittal, "gg")
})

test_that("plot_ortho works with explicit coord", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(1000), c(10, 10, 10)), sp)
  pdf(file = NULL)
  on.exit(dev.off(), add = TRUE)
  result <- plot_ortho(vol, coord = c(5L, 5L, 5L), assemble = FALSE)
  expect_true(is.list(result))
  expect_length(result, 3L)
})

test_that("plot_ortho works with crosshair=FALSE and annotate=FALSE", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(1000), c(10, 10, 10)), sp)
  pdf(file = NULL)
  on.exit(dev.off(), add = TRUE)
  result <- plot_ortho(vol, crosshair = FALSE, annotate = FALSE, assemble = FALSE)
  expect_s3_class(result$axial, "gg")
  expect_null(tile_crosshair(result$axial))
  expect_length(tile_orientation_letters(result$axial), 0L)
})

test_that("plot_ortho works with range='data'", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(1000), c(10, 10, 10)), sp)
  pdf(file = NULL)
  on.exit(dev.off(), add = TRUE)
  result <- plot_ortho(vol, range = "data", assemble = FALSE)
  expect_s3_class(result$axial, "gg")
})

test_that("plot_ortho works with downsample > 1", {
  sp <- NeuroSpace(c(20L, 20L, 20L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(8000), c(20, 20, 20)), sp)
  pdf(file = NULL)
  on.exit(dev.off(), add = TRUE)
  result <- plot_ortho(vol, downsample = 2L, assemble = FALSE)
  expect_true(is.list(result))
  # Decimation by 2 halves each in-plane dimension of the raster.
  expect_equal(nrow(result$axial$data), 10L * 10L)
})

test_that("plot_ortho unit='mm' path runs without error", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(2, 2, 2))
  vol <- DenseNeuroVol(array(rnorm(1000), c(10, 10, 10)), sp)
  pdf(file = NULL)
  on.exit(dev.off(), add = TRUE)
  # coord in mm — coord_to_grid is available in neuroim2
  result <- plot_ortho(vol, coord = c(5, 5, 5), unit = "mm", assemble = FALSE)
  expect_true(is.list(result))
})

test_that("plot_ortho supports draw=FALSE and validates coordinates", {
  sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
  vol <- DenseNeuroVol(array(rnorm(1000), c(10, 10, 10)), sp)

  result <- plot_ortho(vol, coord = c(5L, 5L, 5L), draw = FALSE, style = "dark",
                       assemble = FALSE)
  expect_true(is.list(result))
  expect_s3_class(result$axial, "gg")
  expect_equal(attr(result, "labels")$title, NULL)

  expect_error(
    plot_ortho(vol, coord = c(5L, 5L), draw = FALSE),
    "`coord`"
  )
  expect_error(
    plot_ortho(vol, coord = c(5L, 5L, 99L), draw = FALSE),
    "`coord`"
  )
})
