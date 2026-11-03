library(neuroim2)

# #62: signed overlays honour the host's diverging ramp
# (options(neuro.palette.diverging)) when ov_cmap is NULL.
# #63: the colorbar samples the overlay palette exactly as the voxels do.

host_ramp <- c("#2856cf", "#8fa8e8", "#f6f5f1", "#eaa092", "#d33f29")
signed_values <- c(-12, -6, -1, 0, 2, 5, 12)

# Colour the overlay gives each value (as overlay_slice_grob() does).
overlay_colours <- function(v, sc, thresh) {
  pos <- neuroim2:::overlay_positions(v, sc$lim, thresh = thresh,
                                      diverging = sc$diverging, skip = sc$skip)
  toupper(sc$pal[neuroim2:::palette_index(pos, length(sc$pal))])
}

# The colorbar's bins: the value at each bin's centre and the colour drawn.
bar_bins <- function(sc, thresh) {
  cb <- neuroim2:::neuro_colorbar(sc$lim, sc$pal, thresh = thresh,
                                  diverging = sc$diverging, skip = sc$skip,
                                  tokens = neuroim2:::neuro_style_tokens("report"))
  d <- cb$data
  shown <- d$xmin == 0  # the sub-threshold band is drawn as a narrow strip
  data.frame(value = ((d$ymin + d$ymax) / 2)[shown], fill = toupper(d$fill[shown]))
}

# Colour the colorbar shows at each value.
bar_colours <- function(v, sc, thresh) {
  b <- bar_bins(sc, thresh)
  half <- diff(sc$lim) / (2 * 256)
  vapply(v, function(x) b$fill[which(abs(b$value - x) <= half + 1e-9)[[1L]]], "")
}

# Colours equal up to the rounding of Lab interpolation.
near <- function(a, b) {
  all(abs(grDevices::col2rgb(substr(a, 1, 7)) - grDevices::col2rgb(substr(b, 1, 7))) <= 3)
}

test_that("a signed overlay uses neuro.palette.diverging when ov_cmap is NULL (#62)", {
  withr::local_options(neuro.palette.diverging = host_ramp)
  sc <- neuroim2:::overlay_scale(signed_values, thresh = 3)
  expect_identical(sc$cmap, host_ramp)
  expect_length(sc$pal, 256L)
  expect_true(near(sc$pal[c(1L, 256L)], host_ramp[c(1L, 5L)]))
  # The ends of the scale take the ends of the host ramp.
  expect_true(near(overlay_colours(c(-sc$lim[2], sc$lim[2]), sc, 3),
                   host_ramp[c(1L, 5L)]))
})

test_that("the host ramp starts at its first colour past the neutral middle (#62)", {
  withr::local_options(neuro.palette.diverging = host_ramp)
  sc <- neuroim2:::overlay_scale(signed_values, thresh = 3)
  at_thresh <- overlay_colours(c(-3, 3), sc, 3)
  expect_true(near(at_thresh[1], host_ramp[2]))
  expect_true(near(at_thresh[2], host_ramp[4]))
})

test_that("without the option, or for unsigned data, the defaults are unchanged (#62)", {
  withr::local_options(neuro.palette.diverging = NULL)
  expect_identical(neuroim2:::overlay_scale(signed_values, thresh = 3)$cmap, "cold_hot")
  expect_identical(neuroim2:::overlay_scale(signed_values, thresh = 3)$skip, 0.3)

  withr::local_options(neuro.palette.diverging = host_ramp)
  expect_identical(neuroim2:::overlay_scale(abs(signed_values), thresh = 3)$cmap, "hot")
  expect_identical(neuroim2:::overlay_scale(signed_values, ov_cmap = "RdBu")$cmap, "RdBu")
})

test_that("an unusable neuro.palette.diverging is reported and ignored (#62)", {
  withr::local_options(neuro.palette.diverging = c("not-a-colour", "red", "blue"))
  expect_warning(sc <- neuroim2:::overlay_scale(signed_values, thresh = 3),
                 "neuro.palette.diverging")
  expect_identical(sc$cmap, "cold_hot")
})

test_that("plot_overlay() draws a signed map in the host ramp (#62)", {
  withr::local_options(neuro.palette.diverging = host_ramp)
  bg <- read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
  stat <- bg
  q <- stats::quantile(bg[bg[] > 0], c(0.3, 0.4, 0.9))
  stat[] <- (bg[] > q[3]) * 5 - (bg[] > q[1] & bg[] < q[2]) * 5  # both signs
  tiles <- plot_overlay(bg, stat, zlevels = 24L, ov_thresh = 1, assemble = FALSE)
  layer <- tiles[[1]]$layers[[which(vapply(tiles[[1]]$layers, function(l)
    inherits(l$geom, "GeomCustomAnn"), logical(1)))[[1L]]]]
  r <- layer$geom_params$grob$raster
  drawn <- unique(toupper(as.vector(r)[substr(as.vector(r), 8, 9) != "00"]))
  expect_true(any(vapply(drawn, near, logical(1), b = host_ramp[5])))
  expect_true(any(vapply(drawn, near, logical(1), b = host_ramp[1])))
  expect_false(any(vapply(drawn, near, logical(1), b = "#ffe45c")))  # cold_hot's top
})

# A guard: the bar and the voxels share one palette rule (palette_index()).
test_that("the colorbar shows the overlay's colour at every value (#63)", {
  ramp <- c(grDevices::colorRampPalette(c("#0a2a8f", "#c9d6f5"))(128),
            grDevices::colorRampPalette(c("#f5d0c9", "#8f0a0a"))(128))
  for (cmap in list(ramp, "cold_hot", host_ramp)) {
    sc <- neuroim2:::overlay_scale(c(-20, 20), ov_range = c(-10.2, 10.2),
                                   thresh = 3.59, ov_cmap = cmap)
    b <- bar_bins(sc, 3.59)
    expect_equal(b$fill, overlay_colours(b$value, sc, 3.59))
    # ... including the saturated ends at the cap.
    expect_equal(b$fill[c(1L, nrow(b))], toupper(sc$pal[c(1L, 256L)]))
  }
})

test_that("a palette passed as colours looks as it does by name", {
  # Only a short anchor ramp (the host shape) starts at its first anchor; a
  # full palette vector keeps the named palette's start (review of #62).
  vals <- c(-10, -6, -3.1, 3.1, 6, 10)
  by_name <- neuroim2:::overlay_scale(vals, ov_range = c(-10, 10), thresh = 3,
                                      ov_cmap = "cold_hot")
  by_vec <- neuroim2:::overlay_scale(vals, ov_range = c(-10, 10), thresh = 3,
                                     ov_cmap = resolve_cmap("cold_hot", 256))
  expect_equal(overlay_colours(vals, by_vec, 3), overlay_colours(vals, by_name, 3))
})

test_that("a three-colour ramp stays graded above the threshold", {
  vals <- c(-10, -6, -3.1, 3.1, 6, 10)
  sc <- neuroim2:::overlay_scale(vals, ov_range = c(-10, 10), thresh = 3,
                                 ov_cmap = c("blue", "white", "red"))
  cols <- overlay_colours(vals, sc, 3)
  expect_length(unique(cols[1:3]), 3L)
  expect_length(unique(cols[4:6]), 3L)
  # ... and the two signs stay distinct at the threshold.
  expect_false(near(cols[3], cols[4]))
})

test_that("report montages box their slice labels and leave L/R to the key (#63)", {
  bg <- read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
  stat <- bg
  stat[] <- as.numeric(bg[] > stats::quantile(bg[], 0.9)) * 3
  geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[[1L]], "")

  report <- plot_overlay(bg, stat, zlevels = c(20L, 26L), ov_thresh = 1,
                         style = "report", assemble = FALSE)
  lab <- built_layers_of(report[[1]], "GeomLabel")
  expect_length(lab, 1L)
  expect_match(as.character(lab[[1]]$label), "^z = ")
  # The box takes the tile's own colour, so it reads as part of the tile.
  expect_equal(toupper(lab[[1]]$fill), toupper(resolve_cmap("grays", 2L)[[1L]]))
  expect_length(tile_orientation_letters(report[[1]]), 0L)

  forced <- plot_overlay(bg, stat, zlevels = 20L, ov_thresh = 1, style = "report",
                         annotate = TRUE, assemble = FALSE)
  expect_setequal(unname(tile_orientation_letters(forced[[1]])), c("L", "R"))

  light <- plot_overlay(bg, stat, zlevels = 20L, ov_thresh = 1, assemble = FALSE)
  expect_setequal(unname(tile_orientation_letters(light[[1]])), c("L", "R"))
  no_key <- plot_overlay(bg, stat, zlevels = 20L, ov_thresh = 1, style = "report",
                         legend = FALSE, assemble = FALSE)
  expect_setequal(unname(tile_orientation_letters(no_key[[1]])), c("L", "R"))
})
