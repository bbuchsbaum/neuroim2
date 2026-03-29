## Snapshot tests for neuroim2 plot functions using vdiffr.
## Each test captures the rendered output and compares against a stored SVG.

library(neuroim2)

test_that("plot(vol) produces expected multi-slice output", {
  skip_if_not_installed("vdiffr")
  set.seed(42)
  vol <- make_vol(c(20L, 20L, 20L))
  vdiffr::expect_doppelganger("vol-basic", function() {
    p <- plot(vol)
    print(p)
  })
})

test_that("plot_ortho produces expected three-plane view", {
  skip_if_not_installed("vdiffr")
  set.seed(42)
  vol <- make_vol(c(20L, 20L, 20L))
  vdiffr::expect_doppelganger("ortho-basic", function() {
    plot_ortho(vol)
  })
})

test_that("plot_ortho without crosshairs and annotations", {
  skip_if_not_installed("vdiffr")
  set.seed(42)
  vol <- make_vol(c(20L, 20L, 20L))
  vdiffr::expect_doppelganger("ortho-no-crosshair", function() {
    plot_ortho(vol, crosshair = FALSE, annotate = FALSE)
  })
})

test_that("plot_montage produces expected faceted output", {
  skip_if_not_installed("vdiffr")
  set.seed(42)
  vol <- make_vol(c(20L, 20L, 20L))
  vdiffr::expect_doppelganger("montage-basic", function() {
    p <- plot_montage(vol, zlevels = c(5L, 10L, 15L), ncol = 3L)
    print(p)
  })
})

test_that("plot_overlay with binary alpha mode produces expected output", {
  skip_if_not_installed("vdiffr")
  set.seed(42)
  bg  <- make_vol(c(20L, 20L, 20L))
  ov  <- make_vol(c(20L, 20L, 20L))
  vdiffr::expect_doppelganger("overlay-binary", function() {
    plot_overlay(bg, ov, zlevels = c(5L, 10L, 15L), ncol = 3L)
  })
})

test_that("plot_overlay with proportional alpha mode produces expected output", {
  skip_if_not_installed("vdiffr")
  set.seed(42)
  bg  <- make_vol(c(20L, 20L, 20L))
  ov  <- make_vol(c(20L, 20L, 20L))
  vdiffr::expect_doppelganger("overlay-proportional", function() {
    plot_overlay(bg, ov, zlevels = c(5L, 10L, 15L),
                 ov_alpha_mode = "proportional", ncol = 3L)
  })
})

test_that("plot(NeuroSlice) produces expected 2D output", {
  skip_if_not_installed("vdiffr")
  set.seed(42)
  vol <- make_vol(c(20L, 20L, 20L))
  sl  <- slice(vol, 10L, along = 3L)
  vdiffr::expect_doppelganger("neuroslice-basic", function() {
    p <- plot(sl)
    print(p)
  })
})
