# Slice labels sit in the top-left corner of a tile. A tile is cropped to the
# head with a margin proportional to the head, while a label has a fixed
# physical size, so on small tiles the label used to land on the anatomy.
# The layout now raises the top of the tile by just enough to clear it.

headroom_anat <- function() {
  read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
}

test_that("no headroom is added when the crop margin already fits the label", {
  # 100 x 100 window drawn 3 in wide: 0.03 in per unit. A 3.1 mm label on a
  # box is about 0.19 in with its inset, i.e. under 10 units.
  expect_equal(label_headroom(100, 100, room = 20, label_size = 3.1,
                              w_max = 3, h_max = 3), 0)
})

test_that("headroom makes exactly enough room when the tile is width-limited", {
  xw <- 100; yw <- 100; room <- 1; size <- 3.1
  e <- label_headroom(xw, yw, room, size, w_max = 1, h_max = 10)
  expect_gt(e, 0)
  ipu <- 1 / xw
  label_in <- (size * ggplot2::.pt + 3.2 + 1.5) / 72
  # distance from foreground to window top, in inches, equals inset + label
  expect_equal((room + e) * ipu, 0.035 * 100 * ipu + label_in)
})

test_that("headroom is solved consistently when the tile is height-limited", {
  xw <- 100; yw <- 100; room <- 1; size <- 3.1
  e <- label_headroom(xw, yw, room, size, w_max = 10, h_max = 1)
  expect_gt(e, 0)
  ipu <- 1 / (yw + e)            # the raised tile exactly fills h_max
  label_in <- (size * ggplot2::.pt + 3.2 + 1.5) / 72
  expect_equal((room + e) * ipu, 0.035 * 100 * ipu + label_in)
})

test_that("smaller tiles and larger labels need more headroom", {
  small <- label_headroom(100, 100, 2, 3.1, w_max = 1, h_max = 10)
  large <- label_headroom(100, 100, 2, 3.1, w_max = 2, h_max = 10)
  bigger_text <- label_headroom(100, 100, 2, 5, w_max = 1, h_max = 10)
  expect_gt(small, large)
  expect_gt(bigger_text, small)
  # an unboxed label (the montage) is a little shorter
  expect_lt(label_headroom(100, 100, 2, 3.1, w_max = 1, h_max = 10, boxed = FALSE), small)
})

test_that("degenerate inputs add nothing", {
  expect_equal(label_headroom(100, 100, NA_real_, 3.1, 1, 1), 0)
  expect_equal(label_headroom(100, 100, NULL, 3.1, 1, 1), 0)
  expect_equal(label_headroom(100, 0, 1, 3.1, 1, 1), 0)
})

test_that("the crop window records the top of the foreground", {
  anat <- headroom_anat()
  win <- foreground_crop_window(anat, zlevels = c(20, 30), along = 3L)
  top <- attr(win, "fg_top")
  expect_true(is.numeric(top) && length(top) == 1L)
  expect_lte(top, win$ylim[2] + 1e-8)
  expect_gt(top, mean(win$ylim))
})

test_that("raise_tile_top() lifts the window and the label, not the original", {
  anat <- headroom_anat()
  tiles <- plot_overlay(anat, anat, zlevels = c(20, 30), ov_thresh = 5000,
                        assemble = FALSE)
  p <- tiles[[1]]
  label_y <- function(q) {
    l <- Filter(function(l) inherits(l$geom, "GeomLabel"), q$layers)[[1]]
    l$data$y
  }
  lim0 <- p$coordinates$limits
  y0 <- label_y(p)

  q <- raise_tile_top(p, 7)
  expect_equal(q$coordinates$limits$y, lim0$y + c(0, 7))
  expect_equal(q$coordinates$limits$x, lim0$x)
  expect_equal(label_y(q), y0 + 7)
  expect_no_error(ggplot2::ggplot_build(q))

  # tiles are reused when a figure is re-fitted to another canvas
  expect_equal(p$coordinates$limits, lim0)
  expect_equal(label_y(p), y0)
  expect_identical(raise_tile_top(p, 0), p)
})

test_that("a small canvas raises the tiles and a large one leaves them alone", {
  anat <- headroom_anat()
  tiles <- plot_overlay(anat, anat, zlevels = c(20, 24, 28, 32, 36, 40),
                        ov_thresh = 5000, assemble = FALSE)
  win <- foreground_crop_window(anat, c(20, 24, 28, 32, 36, 40), 3L,
                                extra = list(anat), extra_thresh = 5000)
  room <- list(room = win$ylim[2] - attr(win, "fg_top"), label_size = 3.1)
  aspect <- diff(win$xlim) / diff(win$ylim)
  tok <- neuro_style_tokens("light")
  # highest tile top anywhere in the assembled figure
  top_of <- function(fig) {
    tops <- numeric()
    walk <- function(x) {
      y <- tryCatch(x$coordinates$limits$y, error = function(e) NULL)
      if (is.numeric(y) && length(y) == 2L) tops <<- c(tops, y[2L])
      kids <- tryCatch(x$patches$plots, error = function(e) NULL)
      for (k in kids) walk(k)
    }
    walk(fig)
    max(tops)
  }
  small <- neuro_assemble(tiles, ncol = 6, tile_aspect = aspect, tokens = tok,
                          canvas = c(7.5, 2), label_room = room)
  off <- neuro_assemble(tiles, ncol = 6, tile_aspect = aspect, tokens = tok,
                        canvas = c(7.5, 2), label_room = NULL)
  expect_gt(top_of(small), top_of(off))
  expect_equal(top_of(off), win$ylim[2])
})

test_that("plot_montage() and plot_overlay() still build at README size", {
  anat <- headroom_anat()
  expect_no_error(plot_montage(anat, zlevels = seq(-30, 60, by = 18), unit = "mm",
                               ncol = 6, canvas = c(7.5, 1.9)))
  expect_no_error(plot_overlay(anat, anat, zlevels = c(40, 52, 64), unit = "mm",
                               ov_thresh = 5000, canvas = c(7.5, 2.1)))
})
