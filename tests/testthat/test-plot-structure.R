## Structural checks for the plotting entry points.
##
## These replace the vdiffr golden-image comparisons that used to live in
## test-plot-snapshots.R. Those compared SVG text byte-for-byte, which encodes
## the font metrics of the machine that produced the snapshot, so they failed on
## every platform in the R-CMD-check matrix except the one that generated them.
## The harness is preserved under dev/visual-snapshots/ for reviewing
## intentional visual changes on a reference machine.
##
## What is checked here is what a rendering regression actually breaks: that the
## plot builds, that it has the expected number of panels and layers, and that
## the data reaching the renderer has the right shape. None of that depends on
## fonts or on the graphics device.

library(neuroim2)

demo_volume <- function() {
  set.seed(42)
  make_vol(c(20L, 20L, 20L))
}

built_panels <- function(p) {
  b <- ggplot2::ggplot_build(p)
  length(unique(b$data[[1]]$PANEL))
}

test_that("plot(vol) builds a multi-slice montage", {
  p <- plot(demo_volume())

  expect_s3_class(p, "ggplot")
  expect_silent(b <- ggplot2::ggplot_build(p))
  # The default is a 3 x 3 grid of evenly spaced axial slices.
  expect_equal(built_panels(p), 9L)
  expect_gt(nrow(ggplot2::ggplot_build(p)$data[[1]]), 0)
})

test_that("plot(vol) honours the slicing axis", {
  vol <- make_vol(c(8L, 10L, 12L))
  p <- plot(vol, zlevels = c(2L, 6L), along = 1L)

  expect_equal(built_panels(p), 2L)
  # Sagittal slices contain the other two spatial dimensions.
  expect_equal(nrow(ggplot2::ggplot_build(p)$data[[1]]), 2L * 10L * 12L)
  expect_error(plot(vol, along = 4L), "along")
})

test_that("plot(NeuroSlice) builds a single 2D panel", {
  sl <- slice(demo_volume(), 10L, along = 3L)
  p <- plot(sl)

  expect_s3_class(p, "ggplot")
  expect_equal(built_panels(p), 1L)
  # One raster cell per voxel in the slice.
  expect_equal(nrow(ggplot2::ggplot_build(p)$data[[1]]), 20L * 20L)
})

test_that("plot_montage honours zlevels and ncol", {
  p <- plot_montage(demo_volume(), zlevels = c(5L, 10L, 15L), ncol = 3L)

  expect_s3_class(p, "ggplot")
  expect_equal(built_panels(p), 3L)
})

test_that("plot_ortho returns the three anatomical planes", {
  panels <- plot_ortho(demo_volume(), draw = FALSE)

  expect_named(panels, c("axial", "coronal", "sagittal"), ignore.order = TRUE)
  for (p in panels) {
    expect_s3_class(p, "ggplot")
    expect_silent(ggplot2::ggplot_build(p))
  }
})

test_that("plot_ortho drops crosshair and annotation layers when asked", {
  vol <- demo_volume()

  with_extras <- plot_ortho(vol, draw = FALSE)
  without <- plot_ortho(vol, crosshair = FALSE, annotate = FALSE, draw = FALSE)

  # Turning both off must remove layers, not merely restyle them.
  expect_lt(
    length(without$axial$layers),
    length(with_extras$axial$layers)
  )
  expect_silent(ggplot2::ggplot_build(without$axial))
})

test_that("plot_overlay builds in both alpha modes", {
  set.seed(42)
  bg <- make_vol(c(20L, 20L, 20L))
  ov <- make_vol(c(20L, 20L, 20L))

  binary <- plot_overlay(bg, ov, zlevels = c(5L, 10L, 15L), ncol = 3L)
  proportional <- plot_overlay(bg, ov,
    zlevels = c(5L, 10L, 15L),
    ov_alpha_mode = "proportional", ncol = 3L
  )

  for (p in list(binary, proportional)) {
    expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
    expect_silent(print(p))
  }
})

test_that("plot_overlay panels cover the requested slices", {
  set.seed(42)
  bg <- make_vol(c(20L, 20L, 20L))
  ov <- make_vol(c(20L, 20L, 20L))

  panels <- plot_overlay(bg, ov,
    zlevels = c(5L, 10L, 15L), ncol = 3L,
    assemble = FALSE, draw = FALSE
  )

  expect_length(panels, 3L)
  for (p in panels) {
    built <- ggplot2::ggplot_build(p)
    # One raster cell per voxel of the slice, background layer first.
    expect_equal(nrow(built$data[[1]]), 20L * 20L)
  }
})

test_that("plot_overlay rejects an unknown alpha mode", {
  set.seed(42)
  bg <- make_vol(c(20L, 20L, 20L))
  ov <- make_vol(c(20L, 20L, 20L))

  expect_error(
    plot_overlay(bg, ov, zlevels = 10L, ov_alpha_mode = "nonsense", draw = FALSE)
  )
})
