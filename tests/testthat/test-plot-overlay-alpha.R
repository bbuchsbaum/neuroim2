library(testthat)

# Avoid explicit library(neuroim2) to prevent namespace unload conflicts in checks.

test_that("matrix_to_raster_grob accepts per-pixel alpha maps", {
  mat <- matrix(c(0, 1, 2, 3), nrow = 2, ncol = 2)
  amap <- matrix(c(0, 0.25, 0.75, 1), nrow = 2, ncol = 2)

  g <- neuroim2:::matrix_to_raster_grob(
    mat = mat,
    cmap = "grays",
    limits = c(0, 3),
    alpha = 1,
    alpha_map = amap
  )

  expect_s3_class(g, "rastergrob")
  ras <- g[["raster"]]
  expect_s3_class(ras, "raster")
  expect_equal(dim(ras), c(2, 2))

  cols <- as.character(ras)
  expect_true(all(nchar(cols) == 9))

  # Decode alpha from #RRGGBBAA and compare up to raster quantization.
  alpha_vals <- vapply(
    cols,
    function(x) strtoi(substr(x, 8, 9), base = 16L) / 255,
    numeric(1)
  )
  expect_equal(
    unname(sort(alpha_vals)),
    unname(sort(as.numeric(amap))),
    tolerance = 1 / 255
  )
})

test_that("orient_slice_for_raster places voxel (1,1) at bottom-left for identity transform", {
  sp  <- neuroim2::NeuroSpace(c(4, 5, 2))
  vol <- neuroim2::NeuroVol(array(0, dim = c(4, 5, 2)), sp)
  sl  <- neuroim2::slice(vol, 1, along = 3)

  # Marker at voxel (1,1): world coords (0.5, 0.5) → lower-left of ggplot panel
  mat <- matrix(0, nrow = 4, ncol = 5)
  mat[1, 1] <- 1
  # Second marker at voxel (4, 5): upper-right
  mat[4, 5] <- 2

  out <- neuroim2:::orient_slice_for_raster(sl, mat)
  # rasterGrob renders row 1 at TOP, col 1 at LEFT.
  # voxel (1,1) is lower-left, so it should be at out$mat[nrow, 1]
  expect_equal(out$mat[nrow(out$mat), 1], 1)
  # voxel (4,5) is upper-right, so out$mat[1, ncol]
  expect_equal(out$mat[1, ncol(out$mat)], 2)
  # Extent spans pixel edges: (0, dim) since spacing=1 and origin=0
  expect_equal(out$xmin, 0)
  expect_equal(out$xmax, 4)
  expect_equal(out$ymin, 0)
  expect_equal(out$ymax, 5)
})

test_that("orient_slice_for_raster handles flipped x-axis (radiological LPI-like)", {
  # Build a space where voxel i increases as world x DECREASES (x flipped)
  trans <- diag(c(-1, 1, 1, 1))
  trans[1, 4] <- 3     # so voxel (1,1) -> world (3 - 0.5, 0.5) = (2.5, 0.5)
  sp  <- neuroim2::NeuroSpace(c(4, 5, 2), trans = trans)
  vol <- neuroim2::NeuroVol(array(0, dim = c(4, 5, 2)), sp)
  sl  <- neuroim2::slice(vol, 1, along = 3)

  mat <- matrix(0, nrow = 4, ncol = 5)
  mat[1, 1] <- 1   # voxel (1,1) -> world (2.5, 0.5) = RIGHT-bottom due to flip
  mat[4, 5] <- 2   # voxel (4,5) -> world (-0.5, 4.5) = LEFT-top

  out <- neuroim2:::orient_slice_for_raster(sl, mat)
  # voxel (1,1) is bottom-right -> oriented matrix [nrow, ncol]
  expect_equal(out$mat[nrow(out$mat), ncol(out$mat)], 1)
  # voxel (4,5) is top-left -> oriented matrix [1, 1]
  expect_equal(out$mat[1, 1], 2)
})

test_that("plot_overlay tolerates NAs in overlay with ov_thresh > 0", {
  sp <- neuroim2::NeuroSpace(c(5, 5, 3))
  bg <- neuroim2::NeuroVol(array(0, dim = c(5, 5, 3)), sp)
  ov_arr <- array(stats::rnorm(5 * 5 * 3), dim = c(5, 5, 3))
  ov_arr[1, 1, 1] <- NA_real_
  ov_arr[2, 2, 2] <- NA_real_
  ov <- neuroim2::NeuroVol(ov_arr, sp)

  tf <- tempfile(fileext = ".png")
  grDevices::png(tf)
  on.exit({ grDevices::dev.off(); unlink(tf) }, add = TRUE)

  expect_silent(
    neuroim2::plot_overlay(
      bgvol = bg, overlay = ov, zlevels = c(1, 2), along = 3,
      ov_thresh = 0.5, ncol = 2
    )
  )
  expect_silent(
    neuroim2::plot_overlay(
      bgvol = bg, overlay = ov, zlevels = c(1, 2), along = 3,
      ov_alpha_mode = "proportional", ov_thresh = 0.5, ncol = 2
    )
  )
})

test_that("plot_overlay supports ov_alpha_mode='proportional'", {
  sp <- neuroim2::NeuroSpace(c(6, 6, 4))
  bg <- neuroim2::NeuroVol(array(0, dim = c(6, 6, 4)), sp)
  ov <- neuroim2::NeuroVol(array(stats::rnorm(6 * 6 * 4), dim = c(6, 6, 4)), sp)

  tf <- tempfile(fileext = ".png")
  grDevices::png(tf)
  on.exit({
    grDevices::dev.off()
    unlink(tf)
  }, add = TRUE)

  expect_silent(
    neuroim2::plot_overlay(
      bgvol = bg,
      overlay = ov,
      zlevels = c(1, 2, 3),
      along = 3,
      ov_alpha_mode = "proportional",
      ov_thresh = 0,
      ncol = 3
    )
  )
})
