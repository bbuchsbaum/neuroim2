#' Composite a statistical map on a structural background
#'
#' Draws a grid of slices through a structural background (e.g. a T1) with a
#' thresholded statistical map on top, in the style of a journal figure: black
#' tiles cropped to the head, world-coordinate slice labels, L/R markers, and a
#' compact colorbar that marks the threshold.
#'
#' @param bgvol Background 3D volume.
#' @param overlay Overlay 3D volume on the same NeuroSpace grid as `bgvol`.
#' @param zlevels Slices to plot, as indices along `along` (\code{unit =
#'   "index"}, the default) or world coordinates in mm (\code{unit = "mm"}).
#'   \code{NULL} (default) chooses \code{n_slices} slices spread over the
#'   extent of the supra-threshold overlay.
#' @param along Native voxel-grid axis for slicing. Display orientation and
#'   anatomical plane labels are inferred from the image affine.
#' @param bg_cmap Background palette (e.g., "grays").
#' @param ov_cmap Overlay palette. \code{NULL} (default) chooses automatically:
#'   the two-sided \code{"cold_hot"} map for signed data, \code{"hot"}
#'   otherwise. Any name accepted by [resolve_cmap()] or a vector of colours.
#' @param bg_range,ov_range Background/overlay scaling. Either a mode string,
#'   \code{"robust"} or \code{"data"}, or an explicit numeric \code{c(lo, hi)}
#'   to pin the scale (e.g. \code{ov_range = c(-6, 6)}) for consistent
#'   colouring across figures and subjects. The robust background window is
#'   computed over head voxels only; the robust overlay scale is computed once
#'   from the whole overlay volume (supra-threshold values when a threshold is
#'   set) and rounded to two significant digits, so a given map always gets
#'   the same scale whichever slices are shown.
#' @param probs Quantiles for robust scaling.
#' @param ov_thresh Numeric threshold; values with |v| < thresh are not drawn.
#' @param ov_alpha Global opacity of the overlay (0..1).
#' @param ov_alpha_mode One of \code{"binary"} (default: every supra-threshold
#'   voxel fully opaque), \code{"proportional"}, \code{"ramp"}, or
#'   \code{"soft"} (opacity rises nonlinearly with magnitude; the curve
#'   \code{alpha = floor + (1 - floor) * t^gamma} self-tunes \code{gamma} from
#'   the data). In the graded modes every voxel that passes a threshold keeps at
#'   least 60\% opacity, and the colorbar is faded with the same curve so it
#'   matches the picture.
#' @param alpha_gamma Optional exponent for \code{ov_alpha_mode = "soft"}.
#'   \code{NULL} (default) auto-tunes it from the data.
#' @param ov_symmetric Logical or \code{NULL}. \code{NULL} (default) uses
#'   symmetric limits around zero when the overlay has both signs.
#' @param ov_cap Optional numeric; the magnitude at the upper end of the
#'   colour/opacity scale. Defaults to the data-driven limit.
#' @param ncol Number of columns. \code{NULL} (default) picks the layout that
#'   best fills the canvas.
#' @param title,subtitle,caption Optional figure labels, left-aligned with the
#'   tiles.
#' @param draw Logical; if \code{TRUE}, also print the figure immediately (and
#'   return it invisibly). By default the figure is returned visibly, like a
#'   ggplot, so it prints at the console or in a knitr chunk.
#' @param style Visual style: \code{"light"} (white card), \code{"report"}
#'   (warm off-white card with the key strip on), or \code{"dark"} (black
#'   card). Tiles are black in every style.
#' @param enhance Display-only enhancement of the (unsmoothed) statistical
#'   \code{overlay}. \code{FALSE} (default) leaves it untouched; \code{TRUE}
#'   applies \code{\link{enhance_stat_map}} with defaults; a named \code{list}
#'   is forwarded as arguments to \code{enhance_stat_map()}.
#' @param assemble Logical; if \code{TRUE} (default), return one assembled
#'   \pkg{patchwork} figure. If \code{FALSE}, return the list of per-slice
#'   ggplots.
#' @param colorbar Logical; draw the colorbar (default \code{TRUE}).
#' @param cbar_title Character; the quantity label drawn above the colorbar.
#'   Defaults to \code{"value"}; set it to the statistic actually shown (e.g.
#'   \code{"t"}, \code{"Semipartial r"}).
#' @param legend Logical or \code{NULL}; add a one-line key under the tiles
#'   (plane, neurological convention, and what the threshold shows).
#'   \code{NULL} (default) shows it for \code{style = "report"} only.
#' @param crop Logical; crop every panel to the head bounding box (shared
#'   across slices, and always containing every supra-threshold voxel).
#' @param interpolate Logical; smooth the background raster (default
#'   \code{TRUE}). The overlay itself is always drawn with crisp voxels.
#' @param unit \code{"index"} or \code{"mm"}: how \code{zlevels} is interpreted.
#'   Panels are always labelled in world coordinates (mm).
#' @param annotate Logical; draw L/R orientation letters on the first panel.
#' @param n_slices Number of slices chosen when \code{zlevels} is \code{NULL}.
#' @param canvas Optional \code{c(width, height)} in inches to freeze the
#'   layout for one size. By default (\code{NULL}) the layout is re-fitted to
#'   whatever device the figure is drawn on, including \code{ggsave()}.
#'
#' @return A figure (class \code{neuro_fig}, a \pkg{patchwork} whose
#'   layout is re-fitted to the device it is drawn on) when
#'   \code{assemble = TRUE} or a list of
#'   ggplots (\code{assemble = FALSE}); invisibly when \code{draw = TRUE}.
#'
#' @details
#' \strong{Signed maps.} For overlays with both signs (t/z/contrast maps), the
#' default palette is two-sided and the limits symmetric, so negative values
#' are as visible as positive ones. With a threshold, the colour ramp starts at
#' a saturated colour at \eqn{\pm}threshold and brightens toward the cap; the
#' colorbar shows the sub-threshold band in a neutral tone and ticks the
#' threshold and cap. When the data exceed the cap, the end tick reads
#' \eqn{\ge} cap.
#'
#' \strong{Saving.} The figure's layout is fitted to the device it is drawn
#' on: \code{p <- plot_overlay(...); ggsave("fig.png", p, width = 6, height =
#' 9)} re-arranges the tiles for a 6 x 9 in page, and additions such as
#' \code{p + patchwork::plot_annotation(title = "...")} are kept. Pass
#' \code{canvas} only to freeze the layout for one size.
#'
#' @examples
#' \donttest{
#' bg <- read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
#' stat <- bg
#' stat[] <- rnorm(length(stat)) * (bg[] > stats::quantile(bg[], 0.6))
#' p <- plot_overlay(bg, stat, ov_thresh = 1.5, cbar_title = "z")
#' tf <- tempfile(fileext = ".png")
#' ggplot2::ggsave(tf, p, width = 7, height = 5)
#' }
#' @family plot_neuro
#' @export
plot_overlay <- function(
  bgvol, overlay, zlevels = NULL, along = 3L,
  bg_cmap = "grays", ov_cmap = NULL,
  bg_range = c("robust","data"), ov_range = c("robust","data"),
  probs = c(.02,.98), ov_thresh = 0, ov_alpha = 1,
  ov_alpha_mode = c("binary", "proportional", "ramp", "soft"), ov_symmetric = NULL,
  alpha_gamma = NULL,
  ov_cap = NULL, ncol = NULL, title = NULL, subtitle = NULL, caption = NULL,
  draw = FALSE, style = c("light", "dark", "report"), enhance = FALSE,
  assemble = TRUE, colorbar = TRUE, legend = NULL,
  crop = TRUE, interpolate = TRUE, cbar_title = "value",
  unit = c("index", "mm"), annotate = TRUE, n_slices = 12L, canvas = NULL
) {
  unit_missing <- missing(unit)
  assert_same_neuro_grid(bgvol, overlay = overlay)
  cbar_title <- validate_cbar_title(cbar_title)
  ov_alpha_mode <- match_choice(ov_alpha_mode, c("binary", "proportional", "ramp", "soft"),
                                "ov_alpha_mode")
  style <- match_choice(style, c("light", "dark", "report"))
  unit <- match_choice(unit, c("index", "mm"), "unit")
  tokens <- neuro_style_tokens(style)

  is_report   <- identical(style, "report")
  show_legend <- if (is.null(legend)) is_report else isTRUE(legend)
  do_crop     <- if (is.null(crop)) TRUE else isTRUE(crop)
  interp_bg   <- isTRUE(interpolate)

  # Optional display-only enhancement of the (unsmoothed) statistical overlay.
  overlay <- apply_enhance_arg(overlay, enhance)
  ov_thresh <- check_overlay_args(ov_thresh, ov_alpha)
  check_count(n_slices, "n_slices")
  unit_hint(zlevels, unit_missing, bgvol)
  ov_all <- as.numeric(as.array(overlay))

  if (is.null(zlevels)) {
    support <- array(is.finite(ov_all) & ov_all != 0 & abs(ov_all) >= ov_thresh,
                     dim(overlay)[1:3])
    zlevels <- resolve_slice_levels(NULL, bgvol, along, n = n_slices,
                                    support = if (any(support)) support else NULL)
  } else {
    zlevels <- resolve_slice_levels(zlevels, bgvol, along, unit = unit)
  }
  panel_args <- validate_slice_panel_args(zlevels, along, dim(bgvol), 1L)
  zlevels <- panel_args$zlevels
  along <- panel_args$along

  bg_vals <- unlist(lapply(zlevels, function(z) as.numeric(volume_slice_matrix(bgvol, z, along = along))))
  bg_lim <- background_display_limits(bg_range, bg_vals, probs = probs)

  # One colour scale for the whole map (reproducible across slice choices).
  scale <- overlay_scale(ov_all, ov_range = ov_range, probs = probs,
                         thresh = ov_thresh, ov_cmap = ov_cmap,
                         ov_symmetric = ov_symmetric, ov_cap = ov_cap)
  cap <- max(abs(scale$lim))
  soft <- if (ov_alpha_mode == "soft") {
    soft_alpha_params(abs(ov_all), thresh = ov_thresh, cap = cap, gamma = alpha_gamma)
  }
  alpha_fun <- overlay_alpha_fun(ov_alpha_mode, ov_thresh, cap, soft = soft)

  # Shared head bounding box (so every panel is framed identically).
  crop_win <- if (isTRUE(do_crop)) {
    foreground_crop_window(bgvol, zlevels, along, extra = list(overlay),
                           extra_thresh = max(ov_thresh, 0))
  }

  first <- orient_volume_slice_for_raster(bgvol, zlevels[[1L]], along = along)
  plots <- lapply(seq_along(zlevels), function(k) {
    z <- zlevels[[k]]
    bg_oriented <- if (k == 1L) first else orient_volume_slice_for_raster(bgvol, z, along = along)
    g_ov <- overlay_slice_grob(overlay, z, along, scale, ov_thresh, alpha_fun,
                               alpha = ov_alpha)
    neuro_tile(
      bg_oriented, bg_lim = bg_lim, bg_cmap = bg_cmap, layers = list(g_ov),
      window = crop_win, label = slice_world_label(bgvol, z, along),
      orient = if (isTRUE(annotate) && k == 1L) orientation_letters(bg_oriented) else NULL,
      tokens = tokens, interpolate = interp_bg
    )
  })

  cbar <- if (isTRUE(colorbar)) {
    neuro_colorbar(scale$lim, scale$pal, thresh = ov_thresh,
                   diverging = scale$diverging, tokens = tokens, title = cbar_title,
                   alpha_fun = if (ov_alpha_mode == "binary") NULL else alpha_fun,
                   alpha = ov_alpha, over_hi = scale$over_hi, over_lo = scale$over_lo)
  }
  key <- if (isTRUE(show_legend)) {
    overlay_key(ov_thresh, scale$diverging, cbar_title, tokens, plane = first$plane)
  }
  build <- function(cv) neuro_assemble(plots, ncol = ncol, tile_aspect = tile_aspect_of(crop_win, first),
                        colorbar = cbar, key = key, tokens = tokens,
                        title = title, subtitle = subtitle, caption = caption,
                        canvas = cv)
  cv0 <- canvas_size(canvas)
  fig <- build(cv0)
  if (is.null(canvas)) fig <- neuro_figure(fig, build, cv0, tokens)
  neuro_finish(fig, plots, draw = draw, assemble = assemble, title = title,
               subtitle = subtitle, caption = caption, style = style,
               panel_names = vapply(zlevels, function(z) slice_world_label(bgvol, z, along), ""))
}
