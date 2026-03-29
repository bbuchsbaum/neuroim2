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
