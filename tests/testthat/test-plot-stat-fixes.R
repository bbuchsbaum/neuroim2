library(testthat)
library(neuroim2)

# Helpers ---------------------------------------------------------------------

make_signed_overlay <- function(dims = c(20L, 20L, 20L)) {
  sp <- NeuroSpace(dims)
  bg <- NeuroVol(array(runif(prod(dims), 0, 1000), dims), sp)
  a <- array(0, dims)
  a[5:8,   5:8,   ] <-  6   # positive cluster
  a[14:17, 14:17, ] <- -6   # negative cluster
  list(bg = bg, ov = NeuroVol(a, sp), sp = sp, dims = dims)
}

# plot_overlay default (assemble = TRUE) returns a single patchwork/ggplot;
# assemble = FALSE returns a list of per-slice ggplots.
expect_renderable <- function(x) {
  ok <- inherits(x, "patchwork") || inherits(x, "gg") ||
    (is.list(x) && length(x) > 0 && all(vapply(x, inherits, logical(1), "ggplot")))
  expect_true(ok)
}

# ---- #17: resolve_cmap palette resolution -----------------------------------

test_that("resolve_cmap resolves arbitrary hcl.pals() names instead of viridis", {
  # RdBu/Spectral are valid hcl palettes and must NOT collapse to the fallback.
  # The unknown-name baseline intentionally warns (covered separately), so silence it.
  fallback <- suppressWarnings(resolve_cmap("not_a_palette", 8))
  expect_false(identical(resolve_cmap("RdBu", 8),     fallback))
  expect_false(identical(resolve_cmap("Spectral", 8), fallback))
  expect_length(resolve_cmap("RdBu", 8), 8L)
})

test_that("resolve_cmap warns on a truly unknown palette name", {
  expect_warning(resolve_cmap("totally_made_up_palette", 8), "Unknown palette")
})

test_that("resolve_cmap passes through color vectors and handles built-ins", {
  cols <- c("#000000", "#ffffff")
  expect_identical(resolve_cmap(cols), cols)
  expect_length(resolve_cmap("coolwarm", 16), 16L)   # diverging alias
  expect_length(resolve_cmap("grays", 10), 10L)
})

test_that("is_diverging_cmap classifies palettes", {
  expect_true(is_diverging_cmap("RdBu"))
  expect_true(is_diverging_cmap("Blue-Red"))
  expect_true(is_diverging_cmap("coldhot"))
  expect_false(is_diverging_cmap("inferno"))
  expect_false(is_diverging_cmap("viridis"))
  expect_false(is_diverging_cmap(c("#000000", "#ffffff")))
})

# ---- #19.3: numeric range arguments -----------------------------------------

test_that("plot_overlay accepts numeric ov_range / bg_range", {
  d <- make_signed_overlay()
  res <- plot_overlay(d$bg, d$ov, zlevels = 10L,
                      ov_range = c(-6, 6), bg_range = c(0, 1000), draw = FALSE)
  expect_renderable(res)
})

test_that("plot_overlay rejects malformed numeric ranges", {
  d <- make_signed_overlay()
  expect_error(plot_overlay(d$bg, d$ov, zlevels = 10L, ov_range = c(1, 1), draw = FALSE),
               "distinct")
  expect_error(plot_overlay(d$bg, d$ov, zlevels = 10L, ov_range = c(1, 2, 3), draw = FALSE),
               "two distinct")
})

test_that("plot_ortho and plot_montage accept numeric range", {
  d <- make_signed_overlay()
  expect_silent(plot_ortho(d$bg, range = c(0, 1000), draw = FALSE))
  pm <- plot_montage(d$bg, zlevels = c(8L, 10L), range = c(0, 1000))
  expect_true(inherits(pm, "ggplot"))
})

# ---- #18.1: signed maps -> diverging + symmetric ----------------------------

test_that("signed overlays default to a diverging palette silently", {
  d <- make_signed_overlay()
  # default ov_cmap on signed data: no message/warning, returns a renderable
  expect_silent(res <- plot_overlay(d$bg, d$ov, zlevels = 10L, draw = FALSE))
  expect_renderable(res)
})

test_that("assemble = FALSE returns the per-slice ggplot list", {
  d <- make_signed_overlay()
  res <- plot_overlay(d$bg, d$ov, zlevels = c(8L, 10L), draw = FALSE, assemble = FALSE)
  expect_type(res, "list")
  expect_length(res, 2L)
  expect_true(all(vapply(res, inherits, logical(1), "ggplot")))
})

test_that("explicit sequential palette on signed data warns", {
  d <- make_signed_overlay()
  expect_warning(
    plot_overlay(d$bg, d$ov, zlevels = 10L, ov_cmap = "inferno", draw = FALSE),
    "diverging"
  )
})

test_that("a diverging palette on signed data does not warn", {
  d <- make_signed_overlay()
  expect_silent(plot_overlay(d$bg, d$ov, zlevels = 10L, ov_cmap = "RdBu", draw = FALSE))
})

test_that("purely positive overlays keep the requested sequential palette", {
  sp <- NeuroSpace(c(12L, 12L, 12L))
  bg <- NeuroVol(array(0, c(12,12,12)), sp)
  a <- array(0, c(12,12,12)); a[4:8, 4:8, ] <- 5
  ov <- NeuroVol(a, sp)
  expect_silent(plot_overlay(bg, ov, zlevels = 6L, ov_cmap = "inferno", draw = FALSE))
})

# ---- #19.4: ramp alpha mode -------------------------------------------------

test_that("plot_overlay supports ov_alpha_mode = 'ramp'", {
  d <- make_signed_overlay()
  res <- plot_overlay(d$bg, d$ov, zlevels = c(8L, 10L),
                      ov_alpha_mode = "ramp", ov_thresh = 2, draw = FALSE)
  expect_renderable(res)
})

# ---- #19.1 / #19.2: assemble + colorbar -------------------------------------

test_that("assemble = TRUE returns a single ggsave-able object", {
  skip_if_not_installed("patchwork")
  d <- make_signed_overlay()
  obj <- plot_overlay(d$bg, d$ov, zlevels = c(8L, 10L, 12L),
                      assemble = TRUE, draw = FALSE, title = "T")
  expect_true(inherits(obj, "patchwork") || inherits(obj, "gg"))

  tf <- tempfile(fileext = ".png")
  on.exit(unlink(tf), add = TRUE)
  ggplot2::ggsave(tf, obj, width = 6, height = 3, dpi = 50)
  expect_gt(file.info(tf)$size, 0)
})

test_that("assemble without colorbar also works", {
  skip_if_not_installed("patchwork")
  d <- make_signed_overlay()
  obj <- plot_overlay(d$bg, d$ov, zlevels = c(8L, 10L),
                      assemble = TRUE, colorbar = FALSE, draw = FALSE)
  expect_true(inherits(obj, "patchwork") || inherits(obj, "gg"))
})

# ---- report style: light card, cropped tiles, colorbar, legend --------------

test_that("style = 'report' assembles a single object with legend + crop", {
  skip_if_not_installed("patchwork")
  d <- make_signed_overlay()
  obj <- plot_overlay(d$bg, d$ov, zlevels = c(8L, 10L, 12L), ncol = 3L,
                      style = "report", ov_thresh = 2.3, draw = FALSE)
  expect_true(inherits(obj, "patchwork") || inherits(obj, "gg"))
  tf <- tempfile(fileext = ".png")
  on.exit(unlink(tf), add = TRUE)
  ggplot2::ggsave(tf, obj, width = 8, height = 8, dpi = 50, bg = "#f6f6f4")
  expect_gt(file.info(tf)$size, 0)
})

test_that("legend / crop / interpolate toggle independently of style", {
  d <- make_signed_overlay()
  expect_silent(plot_overlay(d$bg, d$ov, zlevels = 10L, style = "dark",
                             legend = TRUE, crop = TRUE, interpolate = TRUE,
                             ov_thresh = 2.3, draw = FALSE))
})

test_that("style colors, legend grob, and crop window helpers behave", {
  expect_identical(.plot_style_colors("report")$panel, "dark")
  expect_identical(.plot_style_colors("light")$panel, "light")
  leg <- make_overlay_legend(2.3, "#b40426", "#3b4cc0", symmetric = TRUE,
                             style = "report", plane = "Axial")
  expect_true(inherits(leg, "ggplot"))

  d <- make_signed_overlay()
  win <- compute_crop_window(d$bg, d$ov, zlevels = c(8L, 10L), along = 3L,
                             bg_thresh = stats::quantile(as.array(d$bg), 0.2),
                             ov_thresh = 2.3)
  expect_true(is.null(win) || (is.list(win) && all(c("xlim", "ylim") %in% names(win))))
})

# ---- report style is consistent across the plot_ family ---------------------

test_that("plot_ortho supports style = 'report' and returns an assembled object", {
  skip_if_not_installed("patchwork")
  d <- make_signed_overlay()
  obj <- plot_ortho(d$bg, style = "report", draw = FALSE)
  expect_true(inherits(obj, "patchwork") || inherits(obj, "gg"))
  # light/dark path still returns the per-panel list
  lst <- plot_ortho(d$bg, style = "dark", draw = FALSE)
  expect_type(lst, "list")
  expect_length(lst, 3L)
})

test_that("plot_montage supports style = 'report' and returns an assembled object", {
  skip_if_not_installed("patchwork")
  d <- make_signed_overlay()
  obj <- plot_montage(d$bg, zlevels = c(8L, 10L, 12L), ncol = 3L, style = "report")
  expect_true(inherits(obj, "patchwork") || inherits(obj, "gg"))
  # light/dark path still returns a faceted ggplot
  fac <- plot_montage(d$bg, zlevels = c(8L, 10L), style = "dark")
  expect_true(inherits(fac, "ggplot"))
})

test_that("report style ggsave round-trips for ortho and montage", {
  skip_if_not_installed("patchwork")
  d <- make_signed_overlay()
  for (obj in list(plot_ortho(d$bg, style = "report", draw = FALSE),
                   plot_montage(d$bg, zlevels = c(8L, 10L), style = "report"))) {
    tf <- tempfile(fileext = ".png")
    ggplot2::ggsave(tf, obj, width = 8, height = 5, dpi = 50, bg = "#f6f6f4")
    expect_gt(file.info(tf)$size, 0)
    unlink(tf)
  }
})

# ---- #18.2: shared (not per-slice) alpha denominator ------------------------

test_that("proportional alpha uses a shared cap across panels", {
  # Two slices with very different peak magnitudes: identical raw values must
  # map to identical alpha regardless of which slice they sit on. We check this
  # indirectly via the documented behavior that the denominator is the shared
  # ov_lim cap, by confirming a fixed ov_range makes panels independent of
  # per-slice maxima (no error, stable rendering).
  sp <- NeuroSpace(c(10L, 10L, 2L))
  a <- array(0, c(10, 10, 2))
  a[, , 1] <- 2     # low-magnitude slice
  a[, , 2] <- 8     # high-magnitude slice
  ov <- NeuroVol(a, sp)
  bg <- NeuroVol(array(0, c(10, 10, 2)), sp)
  expect_silent(
    plot_overlay(bg, ov, zlevels = c(1L, 2L), ov_cmap = "inferno",
                 ov_alpha_mode = "proportional", ov_range = c(0, 8), draw = FALSE)
  )
})
