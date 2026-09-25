library(testthat)

# Avoid explicit library(neuroim2) to prevent namespace unload conflicts in checks.

# grob_labels() lives in helper-plot.R.

signed_overlay_fixture <- function() {
  sp <- neuroim2::NeuroSpace(c(8L, 8L, 4L), spacing = c(2, 2, 2))
  set.seed(42)
  bg <- neuroim2::NeuroVol(array(runif(8 * 8 * 4), c(8, 8, 4)), sp)
  arr <- array(0, c(8, 8, 4))
  # Both signs must sit on the plotted slice: symmetry is detected from the
  # values actually selected by `zlevels`, not from the whole volume.
  arr[3:5, 3:5, 2] <- 4
  arr[6:7, 6:7, 2] <- -4
  list(bg = bg, overlay = neuroim2::NeuroVol(arr, sp))
}

test_that("plot_overlay forwards cbar_title to the assembled colorbar", {
  fx <- signed_overlay_fixture()
  title <- "Delay coefficient (% signal change)"

  p <- neuroim2::plot_overlay(
    fx$bg, fx$overlay, zlevels = 2L, ov_thresh = 2,
    draw = FALSE, colorbar = TRUE, cbar_title = title
  )

  expect_true(title %in% grob_labels(p))
})

test_that("plot_overlay defaults the colorbar title to 'value'", {
  fx <- signed_overlay_fixture()

  p <- neuroim2::plot_overlay(
    fx$bg, fx$overlay, zlevels = 2L, ov_thresh = 2,
    draw = FALSE, colorbar = TRUE
  )

  expect_true("value" %in% grob_labels(p))
})

test_that("the report key does not assert 'activation' and states the signed threshold", {
  fx <- signed_overlay_fixture()

  # style = "report" is the path that enables the key under the panels.
  p <- neuroim2::plot_overlay(
    fx$bg, fx$overlay, zlevels = 2L, ov_thresh = 2,
    draw = FALSE, style = "report", cbar_title = "Semipartial r"
  )
  labels <- grob_labels(p)

  # An overlay may be a correlation or a coefficient, not a BOLD activation.
  expect_false(any(grepl("activation", labels, ignore.case = TRUE)))
  # The key names the displayed quantity and states the two-sided threshold
  # as a magnitude, so it covers both signs.
  expect_true(any(grepl('group("|", "Semipartial r", "|") ~ phantom() >= "2"', labels, fixed = TRUE)))
  expect_true("Semipartial r" %in% labels)
})

test_that("unsigned overlays get a one-sided suprathreshold key", {
  sp <- neuroim2::NeuroSpace(c(8L, 8L, 4L), spacing = c(2, 2, 2))
  set.seed(7)
  bg <- neuroim2::NeuroVol(array(runif(8 * 8 * 4), c(8, 8, 4)), sp)
  arr <- array(0, c(8, 8, 4))
  arr[3:5, 3:5, 2] <- 4
  ov <- neuroim2::NeuroVol(arr, sp)

  p <- neuroim2::plot_overlay(
    bg, ov, zlevels = 2L, ov_thresh = 2, ov_symmetric = FALSE,
    draw = FALSE, style = "report"
  )
  labels <- grob_labels(p)

  expect_false(any(grepl("activation", labels, ignore.case = TRUE)))
  expect_true(any(grepl('"value" >= "2"', labels, fixed = TRUE)))
  # A one-sided map states a signed (not magnitude) threshold.
  expect_false(any(grepl("|value|", labels, fixed = TRUE)))
})

test_that("plot_overlay rejects a malformed cbar_title", {
  fx <- signed_overlay_fixture()
  args <- list(fx$bg, fx$overlay, zlevels = 2L, ov_thresh = 2, draw = FALSE)

  for (bad in list(NULL, 42, c("a", "b"), NA_character_)) {
    expect_error(
      do.call(neuroim2::plot_overlay, c(args, list(cbar_title = bad))),
      "single non-NA character string"
    )
  }
})

test_that("plot_montage forwards cbar_title to its colorbar", {
  sp <- neuroim2::NeuroSpace(c(8L, 8L, 4L), spacing = c(2, 2, 2))
  set.seed(11)
  vol <- neuroim2::NeuroVol(array(runif(8 * 8 * 4), c(8, 8, 4)), sp)

  # plot_montage() accepts any NeuroVol, a statistic map included, so it needs
  # the same escape hatch as plot_overlay(). Grayscale montages hide the
  # colorbar by default, so request it explicitly.
  p <- neuroim2::plot_montage(vol, zlevels = 2L, style = "report",
                              colorbar = TRUE, cbar_title = "Semipartial r")
  expect_true("Semipartial r" %in% grob_labels(p))

  # A non-grayscale palette shows the colorbar (and its title) by default.
  p <- neuroim2::plot_montage(vol, zlevels = 2L, style = "report",
                              cmap = "viridis", cbar_title = "Semipartial r")
  expect_true("Semipartial r" %in% grob_labels(p))
})

test_that("plot_ortho forwards cbar_title in report style", {
  sp <- neuroim2::NeuroSpace(c(8L, 8L, 8L), spacing = c(2, 2, 2))
  set.seed(13)
  vol <- neuroim2::NeuroVol(array(runif(8 * 8 * 8), c(8, 8, 8)), sp)

  p <- neuroim2::plot_ortho(vol, coord = c(4L, 4L, 4L), draw = FALSE,
                            style = "report", colorbar = TRUE,
                            cbar_title = "Delay coefficient")
  expect_true("Delay coefficient" %in% grob_labels(p))

  # With an overlay the colorbar is informative and shown by default.
  ov <- neuroim2::NeuroVol(array(rnorm(8 * 8 * 8), c(8, 8, 8)), sp)
  p <- neuroim2::plot_ortho(vol, coord = c(4L, 4L, 4L), draw = FALSE,
                            style = "report", overlay = ov,
                            cbar_title = "Delay coefficient")
  expect_true("Delay coefficient" %in% grob_labels(p))
})

test_that("colorbar entry points default to the historical 'value' title", {
  sp <- neuroim2::NeuroSpace(c(8L, 8L, 8L), spacing = c(2, 2, 2))
  set.seed(17)
  vol <- neuroim2::NeuroVol(array(runif(8 * 8 * 8), c(8, 8, 8)), sp)

  expect_true("value" %in%
    grob_labels(neuroim2::plot_montage(vol, zlevels = 2L, style = "report",
                                       colorbar = TRUE)))
  expect_true("value" %in%
    grob_labels(neuroim2::plot_ortho(vol, coord = c(4L, 4L, 4L), draw = FALSE,
                                     style = "report", colorbar = TRUE)))
})

test_that("grayscale structural views omit the uninformative colorbar by default", {
  sp <- neuroim2::NeuroSpace(c(8L, 8L, 8L), spacing = c(2, 2, 2))
  set.seed(19)
  vol <- neuroim2::NeuroVol(array(runif(8 * 8 * 8), c(8, 8, 8)), sp)

  expect_false("value" %in%
    grob_labels(neuroim2::plot_montage(vol, zlevels = 2L, style = "report")))
  expect_false("value" %in%
    grob_labels(neuroim2::plot_ortho(vol, coord = c(4L, 4L, 4L), draw = FALSE,
                                     style = "report")))
})
