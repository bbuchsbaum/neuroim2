library(testthat)

# Tests for behaviour introduced by the shared plotting engine
# (R/plot-engine.R): threshold-aware overlay limits, automatic slice
# selection, world-coordinate labels, assembled figures, and keys.

test_that("overlay_display_limits keeps the threshold inside the scale", {
  set.seed(1)
  # Mostly zeros (outside the mask) plus near-zero noise and a handful of
  # supra-threshold voxels: the robust upper quantile falls below the threshold.
  v <- c(rep(0, 5000), rnorm(2000, sd = 0.5), 3.2, 3.4, -3.3)
  lim <- neuroim2:::overlay_display_limits("robust", v, thresh = 3)
  expect_gte(lim[[2L]], 3)
  expect_lte(lim[[1L]], -3)

  # One-sided map: lower end anchored at zero, upper end still reaches thresh.
  pos <- c(rep(0, 5000), abs(rnorm(2000, sd = 0.5)), 3.1, 3.6)
  lim <- neuroim2:::overlay_display_limits("robust", pos, thresh = 3)
  expect_equal(lim[[1L]], 0)
  expect_gte(lim[[2L]], 3)

  # Zeros do not drive the scale: without a threshold the limits come from the
  # non-zero values only, so the robust window is not collapsed toward zero.
  lim0 <- neuroim2:::overlay_display_limits("robust", c(rep(0, 1e4), 5:10))
  expect_gt(lim0[[2L]], 5)

  # An all-zero map yields a valid, non-degenerate scale.
  expect_equal(neuroim2:::overlay_display_limits("robust", rep(0, 10)), c(0, 1))
})

test_that("automatic slice selection skips empty planes", {
  dims <- c(12L, 12L, 30L)
  arr <- array(0, dims)
  arr[3:10, 3:10, 8:22] <- 100 + runif(8 * 8 * 15)
  vol <- neuroim2::NeuroVol(arr, neuroim2::NeuroSpace(dims))

  z <- neuroim2:::default_slice_levels(vol, along = 3L, n = 9L)
  expect_gt(length(z), 1L)
  expect_true(all(z >= 8L & z <= 22L))
  expect_true(all(vapply(z, function(k) any(arr[, , k] > 0), logical(1))))

  # resolve_slice_levels(NULL, ...) is what the plot_* functions call.
  expect_equal(neuroim2:::resolve_slice_levels(NULL, vol, 3L, n = 9L), z)

  # With a support mask (e.g. the supra-threshold overlay), slices are drawn
  # from its extent rather than from the whole head.
  support <- array(FALSE, dims)
  support[5:7, 5:7, 12:15] <- TRUE
  zs <- neuroim2:::default_slice_levels(vol, along = 3L, n = 4L, support = support)
  expect_true(all(zs >= 12L & zs <= 15L))

  # And the montage built from default slices shows only non-empty planes.
  panels <- neuroim2::plot_montage(vol, n_slices = 6L)
  shown <- unique(panels$data$z)
  expect_true(all(vapply(shown, function(k) any(arr[, , k] > 0), logical(1))))
})

test_that("slice_world_label reports world coordinates in mm", {
  sp <- neuroim2::NeuroSpace(c(10L, 12L, 14L), spacing = c(2, 2, 3),
                             origin = c(-10, -20, -30))
  vol <- neuroim2::NeuroVol(array(0, c(10L, 12L, 14L)), sp)

  expect_equal(neuroim2:::slice_world_label(vol, 5L, along = 3L), "z = -18 mm")
  expect_equal(neuroim2:::slice_world_label(vol, 1L, along = 1L), "x = -10 mm")
  expect_equal(neuroim2:::slice_world_label(vol, 6L, along = 2L), "y = -10 mm")

  # unit = "mm" inverts the labelling: a world position maps back to the slice.
  expect_equal(neuroim2:::resolve_slice_levels(-18, vol, 3L, unit = "mm"), 5L)
  expect_equal(neuroim2:::resolve_slice_levels(c(-10, 2), vol, 1L, unit = "mm"),
               c(1L, 7L))

  # With permuted voxel axes the label names the world axis the slice is
  # normal to, not the native index.
  affine <- diag(4)
  affine[1:3, 1:3] <- matrix(c(0, 1, 0,
                               0, 0, 1,
                               1, 0, 0), nrow = 3L)
  pvol <- neuroim2::NeuroVol(array(0, c(5L, 6L, 7L)),
                             neuroim2::NeuroSpace(c(5L, 6L, 7L), trans = affine))
  expect_equal(neuroim2:::slice_world_label(pvol, 3L, along = 3L), "x = 2 mm")
  expect_equal(neuroim2:::slice_world_label(pvol, 5L, along = 2L), "z = 4 mm")
})

test_that("plot_montage labels panels with world coordinates", {
  sp <- neuroim2::NeuroSpace(c(10L, 12L, 14L), spacing = c(2, 2, 3),
                             origin = c(-10, -20, -30))
  set.seed(3)
  vol <- neuroim2::NeuroVol(array(runif(10 * 12 * 14), c(10L, 12L, 14L)), sp)
  p <- neuroim2::plot_montage(vol, zlevels = c(5L, 9L))
  labs <- unlist(lapply(built_layers_of(p, "GeomText"), function(d) as.character(d$label)))
  expect_true(all(c("z = -18 mm", "z = -6 mm") %in% labs))
})

test_that("plot_ortho(assemble = TRUE) returns one patchwork figure", {
  set.seed(5)
  sp <- neuroim2::NeuroSpace(c(12L, 12L, 12L), spacing = c(2, 2, 2))
  vol <- neuroim2::NeuroVol(array(runif(12^3), c(12L, 12L, 12L)), sp)
  fig <- neuroim2::plot_ortho(vol, coord = c(6L, 6L, 6L), draw = FALSE,
                              title = "Ortho title")
  expect_s3_class(fig, "patchwork")
  labels <- grob_labels(fig)
  expect_true("Ortho title" %in% labels)
  # All three planes are present, each labelled in mm (voxel 6 -> 10 mm).
  for (lab in c("x = 10", "y = 10", "z = 10")) {
    expect_true(any(grepl(paste0("^", lab, "( mm)?$"), labels)), info = lab)
  }

  tf <- tempfile(fileext = ".png")
  on.exit(unlink(tf), add = TRUE)
  ggplot2::ggsave(tf, fig, width = 6, height = 3, dpi = 40)
  expect_gt(file.info(tf)$size, 0)
})

test_that("plot_ortho with an overlay centres on the overlay peak", {
  dims <- c(10L, 11L, 12L)
  sp <- neuroim2::NeuroSpace(dims)
  set.seed(9)
  bg <- neuroim2::NeuroVol(array(runif(prod(dims)), dims), sp)
  ov_arr <- array(0, dims)
  ov_arr[3L, 8L, 9L] <- -7  # peak |value|, negative sign
  ov_arr[6L, 4L, 2L] <- 5
  ov <- neuroim2::NeuroVol(ov_arr, sp)

  panels <- neuroim2::plot_ortho(bg, overlay = ov, draw = FALSE, assemble = FALSE)
  expect_world_label(tile_slice_label(panels$sagittal), "x = 2")
  expect_world_label(tile_slice_label(panels$coronal), "y = 7")
  expect_world_label(tile_slice_label(panels$axial), "z = 8")
  # The overlay is drawn as an extra raster layer over the background.
  expect_true(any(vapply(panels$axial$layers, layer_geom, character(1)) == "GeomCustomAnn"))
})

test_that("registration QC keys name both images", {
  dims <- c(8L, 9L, 5L)
  sp <- neuroim2::NeuroSpace(dims)
  bg <- neuroim2::NeuroVol(array(seq_len(prod(dims)), dims), sp)
  ov <- neuroim2::NeuroVol(array(rev(seq_len(prod(dims))), dims), sp)
  e1 <- neuroim2::NeuroVol(array(seq_len(prod(dims)) %% 3L, dims), sp)
  e2 <- neuroim2::NeuroVol(array(seq_len(prod(dims)) %% 5L, dims), sp)

  edge <- neuroim2::plot_edge_overlay(bg, e1, e2, zlevels = 2L, draw = FALSE,
                                      labels = c("MNI template", "registered EPI"))
  labs <- grob_labels(edge)
  expect_true(all(c("MNI template", "registered EPI") %in% labs))

  edge_nokey <- neuroim2::plot_edge_overlay(bg, e1, e2, zlevels = 2L, draw = FALSE,
                                            labels = c("MNI template", "registered EPI"),
                                            legend = FALSE)
  expect_false(any(c("MNI template", "registered EPI") %in% grob_labels(edge_nokey)))

  chk <- neuroim2::plot_checkerboard(bg, ov, zlevels = 2L, tile = 3L, draw = FALSE,
                                     labels = c("T1", "EPI"))
  key <- grep("^Top-left tile", grob_labels(chk), value = TRUE)
  expect_length(key, 1L)
  expect_match(key, "Top-left tile: T1, alternating with EPI", fixed = TRUE)
  expect_true(any(grepl("^3 mm tiles", grob_labels(chk))))

  chk_nokey <- neuroim2::plot_checkerboard(bg, ov, zlevels = 2L, draw = FALSE,
                                           legend = FALSE)
  expect_false(any(grepl("^Top-left tile", grob_labels(chk_nokey))))

  expect_error(
    neuroim2::plot_edge_overlay(bg, e1, e2, zlevels = 2L, draw = FALSE, labels = "one"),
    "`labels`"
  )
  expect_error(
    neuroim2::plot_checkerboard(bg, ov, zlevels = 2L, draw = FALSE, labels = c("a", "b", "c")),
    "`labels`"
  )
})

test_that("overlay_key states the threshold without 'activation' wording", {
  tokens <- neuroim2:::neuro_style_tokens("report")
  key_text <- function(...) paste(grob_labels(neuroim2:::overlay_key(...)), collapse = " ")

  signed <- key_text(2.5, TRUE, "t-statistic", tokens)
  # Two-sided maps state the threshold as a magnitude (both signs shown).
  # The comparison is plotmath (renders on every device, including pdf()).
  expect_match(signed, 'group("|", "t-statistic", "|") ~ phantom() >= "2.5"', fixed = TRUE)
  expect_match(signed, "Axial slices, neurological view (L = left)", fixed = TRUE)
  expect_false(grepl("activation", signed, ignore.case = TRUE))

  unsigned <- key_text(2.5, FALSE, "value", tokens, plane = "Coronal")
  expect_match(unsigned, '"value" >= "2.5"', fixed = TRUE)
  expect_false(grepl("|", unsigned, fixed = TRUE))
  expect_match(unsigned, "Coronal slices, neurological view (L = left)", fixed = TRUE)

  # No threshold: only the orientation note.
  none <- key_text(0, TRUE, "value", tokens)
  expect_equal(none, "Axial slices, neurological view (L = left)")
})

test_that("device-fitted figures re-fit on draw and keep later additions", {
  dims <- c(12L, 14L, 10L)
  sp <- neuroim2::NeuroSpace(dims)
  set.seed(1)
  bg <- neuroim2::NeuroVol(array(runif(prod(dims), 50, 100), dims), sp)
  ov <- neuroim2::NeuroVol(array(rnorm(prod(dims)), dims), sp)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  p <- neuroim2::plot_overlay(bg, ov, ov_thresh = 1, zlevels = 3:6)
  expect_s3_class(p, "neuro_fig")
  # print() and grid.draw() are registered S3 methods (not just exported).
  expect_true(!is.null(getS3method("print", "neuro_fig", optional = TRUE)))
  expect_true(!is.null(getS3method("grid.draw", "neuro_fig", optional = TRUE,
                                   envir = asNamespace("grid"))))

  # Additions are applied now and replayed when the layout is re-fitted.
  p2 <- p + patchwork::plot_annotation(caption = "added later")
  expect_s3_class(p2, "neuro_fig")
  env <- attr(p2, "neuro_figure")
  refit <- env$build(c(4, 9))
  expect_true("added later" %in% grob_labels(refit))
  expect_true("added later" %in% grob_labels(p2))

  tf <- tempfile(fileext = ".png")
  on.exit(unlink(tf), add = TRUE)
  expect_error(ggplot2::ggsave(tf, p2, width = 4, height = 6, dpi = 40), NA)
  expect_gt(file.info(tf)$size, 0)
})

test_that("slice positions given without a unit trigger a one-time hint on non-trivial affines", {
  dims <- c(10L, 10L, 10L)
  sp <- neuroim2::NeuroSpace(dims, origin = c(-20, -30, -40))
  bg <- neuroim2::NeuroVol(array(runif(prod(dims)), dims), sp)
  old <- options(rlib_message_verbosity = "verbose")
  on.exit(options(old), add = TRUE)
  expect_message(neuroim2::plot_montage(bg, zlevels = c(3L, 5L)), "voxel indices")
  expect_no_message(neuroim2::plot_montage(bg, zlevels = c(-35, -33), unit = "mm"))
})

test_that("plot_checkerboard rejects fractional and non-positive tiles", {
  dims <- c(8L, 8L, 4L)
  sp <- neuroim2::NeuroSpace(dims)
  v <- neuroim2::NeuroVol(array(runif(prod(dims)), dims), sp)
  expect_error(neuroim2::plot_checkerboard(v, v, zlevels = 2L, tile = 2.5), "tile")
  expect_error(neuroim2::plot_checkerboard(v, v, zlevels = 2L, tile = 0), "tile")
})

test_that("plot_ortho views share one height and stay above the field-of-view floor", {
  dims <- c(20L, 24L, 16L)
  sp <- neuroim2::NeuroSpace(dims)
  arr <- array(0, dims)
  arr[4:17, 4:21, 3:14] <- 100   # a head that is longer (A-P) than it is tall (S-I)
  vol <- neuroim2::NeuroVol(arr + runif(prod(dims)), sp)
  panels <- neuroim2::plot_ortho(vol, assemble = FALSE)
  heights <- vapply(panels, function(p) diff(p$coordinates$limits$y), numeric(1))
  expect_equal(unname(heights), rep(heights[[1]], 3), tolerance = 1e-8)
  for (nm in c("sagittal", "coronal")) {
    built <- ggplot2::ggplot_build(panels[[nm]])
    ymin_data <- min(built$data[[1]]$y) - 0.5
    expect_gte(panels[[nm]]$coordinates$limits$y[1], ymin_data - 1e-8)
  }
})

test_that("tile padding uses the background palette's lowest colour", {
  dims <- c(10L, 10L, 6L)
  vol <- neuroim2::NeuroVol(array(runif(600), dims), neuroim2::NeuroSpace(dims))
  p <- neuroim2::plot_ortho(vol, cmap = "viridis", assemble = FALSE)$axial
  expect_identical(tolower(p$theme$panel.background$fill),
                   tolower(neuroim2::resolve_cmap("viridis", 2L)[[1L]]))
})

test_that("plot_edge_overlay validates agree_color and can turn agreement off", {
  dims <- c(8L, 8L, 4L)
  sp <- neuroim2::NeuroSpace(dims)
  bg <- neuroim2::NeuroVol(array(runif(prod(dims)), dims), sp)
  e <- neuroim2::NeuroVol(array(seq_len(prod(dims)) %% 3L, dims), sp)
  expect_error(neuroim2::plot_edge_overlay(bg, e, e, zlevels = 2L, agree_color = "notacolour"),
               "agree_color")
  off <- neuroim2::plot_edge_overlay(bg, e, e, zlevels = 2L, agree_color = NA)
  expect_false("both" %in% grob_labels(off))
  on <- neuroim2::plot_edge_overlay(bg, e, e, zlevels = 2L)
  expect_true("both" %in% grob_labels(on))
})
