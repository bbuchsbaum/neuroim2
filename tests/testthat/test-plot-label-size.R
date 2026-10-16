library(neuroim2)

# Regression tests for #51: the in-tile slice label and orientation letters
# were drawn at a hard-coded ggplot2 size (3.1, letters 3.41) that no argument
# or theme could reach.

label_size_bg <- function() {
  read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
}

# Foreground (non-halo) text sizes of one tile or one faceted montage plot:
# slice labels (multi-character) and orientation letters (one character).
tile_text_sizes <- function(p) {
  d <- do.call(rbind, lapply(tile_text_data(p), function(x) x[, c("label", "size")]))
  lab <- as.character(d$label)
  list(label = unique(d$size[nchar(lab) > 1L]),
       letter = unique(d$size[nchar(lab) == 1L]))
}

expect_tile_sizes <- function(p, size, letters = TRUE) {
  s <- tile_text_sizes(p)
  expect_equal(s$label, size)
  if (letters) expect_equal(s$letter, size * 1.1)
}

test_that("plot_ortho() passes label_size to every view (#51)", {
  bg <- label_size_bg()
  expect_true("label_size" %in% names(formals(plot_ortho)))
  default <- plot_ortho(bg, assemble = FALSE)
  for (p in default) expect_tile_sizes(p, 3.1)
  big <- plot_ortho(bg, label_size = 6, assemble = FALSE)
  expect_named(big, c("axial", "coronal", "sagittal"), ignore.order = TRUE)
  for (p in big) expect_tile_sizes(p, 6)
  # The assembled figure carries the same tiles.
  fig <- plot_ortho(bg, label_size = 6)
  expect_s3_class(fig, "patchwork")
})

test_that("plot_montage() passes label_size to slice labels and letters (#51)", {
  bg <- label_size_bg()
  expect_tile_sizes(plot_montage(bg, n_slices = 4L), 3.1)
  expect_tile_sizes(plot_montage(bg, n_slices = 4L, label_size = 5), 5)
})

test_that("plot_overlay() passes label_size to its tiles (#51)", {
  bg <- label_size_bg()
  stat <- bg
  stat[] <- as.numeric(bg[] > stats::quantile(bg[], 0.9)) * 3
  tiles <- plot_overlay(bg, stat, zlevels = c(20L, 26L), ov_thresh = 1,
                        label_size = 4.5, assemble = FALSE)
  expect_tile_sizes(tiles[[1]], 4.5)            # first tile has L/R letters
  expect_tile_sizes(tiles[[2]], 4.5, letters = FALSE)
})

test_that("registration QC plots pass label_size to their tiles (#51)", {
  bg <- label_size_bg()
  moving <- bg
  moving[] <- bg[c(2:dim(bg)[1], 1), , ]
  cb <- plot_checkerboard(bg, moving, zlevels = c(20L, 26L), label_size = 4,
                          assemble = FALSE)
  expect_tile_sizes(cb[[1]], 4)
  eo <- plot_edge_overlay(bg, bg, moving, zlevels = c(20L, 26L),
                          compute_edges = TRUE, label_size = 4, assemble = FALSE)
  expect_tile_sizes(eo[[1]], 4)
})

test_that("label_size is validated (#51)", {
  bg <- label_size_bg()
  for (bad in list(0, -1, NA_real_, Inf, "big", c(3, 4))) {
    expect_error(plot_ortho(bg, label_size = bad, assemble = FALSE),
                 "label_size", info = deparse(bad))
  }
  expect_error(plot_montage(bg, label_size = 0), "label_size")
  expect_error(plot_overlay(bg, bg, label_size = -2), "label_size")
  expect_error(plot_checkerboard(bg, bg, label_size = NA_real_), "label_size")
  expect_error(plot_edge_overlay(bg, bg, bg, label_size = "x"), "label_size")
})
