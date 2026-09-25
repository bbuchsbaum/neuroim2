#' Checkerboard comparison of two registered volumes
#'
#' Alternates tiles from a background volume and a comparison volume on matched
#' slices. This is useful for visual registration QC.
#'
#' @param bgvol Background/reference 3D volume.
#' @param overlay Comparison 3D volume on the same NeuroSpace grid as `bgvol`.
#' @param zlevels Slices to plot: indices along `along` (\code{unit = "index"})
#'   or world coordinates (\code{unit = "mm"}). \code{NULL} (default) picks
#'   \code{n_slices} slices spread over the brain.
#' @param along Native voxel-grid axis for slicing. Display orientation is
#'   inferred from the image affine.
#' @param tile Tile width in voxels (the key states the equivalent in mm).
#'   \code{NULL} (default) picks a clean size in mm giving about ten tiles
#'   across the head.
#' @param cmap Palette used to render the normalized checkerboard image.
#' @param bg_range,ov_range Intensity scaling of each image: \code{"robust"}
#'   (computed over head voxels), \code{"data"}, or numeric \code{c(lo, hi)}.
#'   Each image is windowed independently so tissue contrast, not a brightness
#'   step, is what differs between tiles.
#' @param probs Quantiles for robust scaling.
#' @param ncol Number of columns (\code{NULL} = chosen to fill the canvas).
#' @param title,subtitle,caption Optional layout-level labels used when drawing.
#' @param draw Logical; if \code{TRUE}, also print the figure immediately (and
#'   return it invisibly). By default the figure is returned visibly.
#' @param style Visual style: \code{"light"}, \code{"dark"} or \code{"report"}.
#' @param labels Length-2 character vector naming the two images in the key.
#' @param legend Logical; draw the key under the panels.
#' @param crop Logical; crop panels to the head bounding box.
#' @param mask_background Logical; show the checker pattern only inside the
#'   head (union of both images' foreground), leaving air black.
#' @param unit \code{"index"} or \code{"mm"}: how \code{zlevels} is
#'   interpreted. Panels are labelled in world coordinates.
#' @param annotate Logical; draw orientation letters on the first panel.
#' @param assemble Logical; return one \pkg{patchwork} figure (default) or the
#'   list of panel ggplots.
#' @param n_slices Number of automatically chosen slices.
#' @param match_intensity Logical; map the comparison image's intensities onto
#'   the reference image's distribution (quantile matching over head voxels)
#'   before interleaving, so a brightness or contrast difference between the
#'   images does not masquerade as misalignment. Default \code{TRUE}.
#' @param canvas Optional \code{c(width, height)} in inches to fit the layout
#'   to (default: the open device).
#' @param focus_brain Logical; interleave tiles only inside the (dilated)
#'   bright-tissue mask of \code{bgvol} and show the fixed image alone on the
#'   scalp, so the seams where registration is judged dominate the figure
#'   instead of the stair-stepped head outline. Default \code{TRUE}.
#' @param interpolate Logical; smooth the rendered checkerboard (default
#'   \code{TRUE}), which also softens the head outline.
#' @return A figure (class \code{neuro_fig}, a \pkg{patchwork} whose
#'   layout is re-fitted to the device it is drawn on) when
#'   \code{assemble = TRUE} or a named list
#'   of panel ggplots (\code{assemble = FALSE}); invisibly when
#'   \code{draw = TRUE}.
#' @examples
#' \donttest{
#' fixed <- read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
#' moving <- fixed
#' moving[] <- sqrt(fixed[c(2:dim(fixed)[1], 1), , ])   # shifted, different contrast
#' p <- plot_checkerboard(fixed, moving)
#' ggplot2::ggsave(tempfile(fileext = ".png"), p, width = 8, height = 6)
#' }
#' @family plot_neuro
#' @export
plot_checkerboard <- function(
    bgvol, overlay, zlevels = NULL, along = 3L, tile = NULL, cmap = "grays",
    bg_range = c("robust", "data"), ov_range = c("robust", "data"),
    probs = c(.02, .98), ncol = NULL,
    title = NULL, subtitle = NULL, caption = NULL, draw = FALSE,
    style = c("light", "dark", "report"), labels = c("fixed", "moving"),
    legend = TRUE, crop = TRUE, mask_background = TRUE,
    unit = c("index", "mm"), annotate = TRUE, assemble = TRUE, n_slices = 6L,
    match_intensity = TRUE, canvas = NULL, interpolate = TRUE, focus_brain = TRUE) {
  unit_missing <- missing(unit)
  assert_same_neuro_grid(bgvol, overlay = overlay)
  style <- match_choice(style, c("light", "dark", "report"))
  unit <- match_choice(unit, c("index", "mm"), "unit")
  tokens <- neuro_style_tokens(style)
  labels <- check_labels2(labels)
  check_count(n_slices, "n_slices")
  unit_hint(zlevels, unit_missing, bgvol)
  zlevels <- resolve_slice_levels(zlevels, bgvol, along, unit = unit, n = n_slices)
  panel_args <- validate_slice_panel_args(zlevels, along, dim(bgvol), 1L)
  zlevels <- panel_args$zlevels
  along <- panel_args$along

  sel <- function(vol) unlist(lapply(zlevels, function(z) as.numeric(volume_slice_matrix(vol, z, along))))
  bg_vals <- sel(bgvol); ov_vals <- sel(overlay)
  bg_lim <- background_display_limits(bg_range, bg_vals, probs = probs)
  ov_lim <- background_display_limits(ov_range, ov_vals, probs = probs)
  bg_fg <- foreground_threshold(bg_vals)
  ov_fg <- foreground_threshold(ov_vals)
  brain_thr <- NA_real_
  if (isTRUE(focus_brain) && is.finite(bg_fg)) {
    brain_thr <- foreground_threshold(bg_vals[is.finite(bg_vals) & bg_vals > bg_fg])
  }
  # Quantile-match the comparison image to the reference over head voxels.
  matcher <- NULL
  if (isTRUE(match_intensity) && is.finite(bg_fg) && is.finite(ov_fg)) {
    ref <- sort(bg_vals[is.finite(bg_vals) & bg_vals > bg_fg])
    src <- ov_vals[is.finite(ov_vals) & ov_vals > ov_fg]
    if (length(ref) > 10L && length(src) > 10L) {
      src_ecdf <- stats::ecdf(src)
      matcher <- function(x) {
        out <- x
        ok <- is.finite(x)
        out[ok] <- stats::quantile(ref, pmin(pmax(src_ecdf(x[ok]), 0), 1), names = FALSE, type = 1)
        out[ok & x <= ov_fg] <- x[ok & x <= ov_fg] * (min(ref) / max(ov_fg, 1e-12))
        out
      }
      ov_lim <- bg_lim
    }
  }

  window <- if (isTRUE(crop)) foreground_crop_window(bgvol, zlevels, along, bg_thresh = bg_fg) else NULL
  if (is.null(tile)) {
    o1 <- orient_volume_slice_for_raster(bgvol, zlevels[[1L]], along = along)
    width_mm <- if (is.null(window)) diff(raster_extent_from_centers(o1$x)) else diff(window$xlim)
    px <- min(spacing(space(bgvol))[-along])
    # About ten tiles across the head, rounded to a clean size in mm.
    target_mm <- width_mm / 10
    nice <- c(4, 5, 6, 8, 10, 12, 15, 16, 20, 25, 30, 40, 50)
    tile_mm <- nice[which.min(abs(nice - target_mm))]
    tile <- max(2L, as.integer(round(tile_mm / px)))
  }
  if (!is.numeric(tile) || length(tile) != 1L || !is.finite(tile) || tile < 1 || tile != round(tile)) {
    cli::cli_abort("{.arg tile} must be a positive whole number of voxels.", call = NULL)
  }
  tile <- as.integer(round(tile))

  plots <- lapply(seq_along(zlevels), function(k) {
    z <- zlevels[[k]]
    bg <- volume_slice_matrix(bgvol, z, along = along)
    ov <- volume_slice_matrix(overlay, z, along = along)
    ov_raw <- ov
    if (!is.null(matcher)) ov[] <- matcher(as.numeric(ov))
    if (!identical(dim(bg), dim(ov))) {
      cli::cli_abort("Sliced volumes must have matching dimensions.", call = NULL)
    }
    bg01 <- matrix(pmax(0, pmin(1, rescale01(as.numeric(bg), bg_lim))), nrow(bg))
    ov01 <- matrix(pmax(0, pmin(1, rescale01(as.numeric(ov), ov_lim))), nrow(ov))
    inside <- (is.finite(bg) & bg > bg_fg) | (is.finite(ov_raw) & ov_raw > ov_fg)
    # Dark CSF between brain and skull falls below the head/air threshold;
    # fill interior holes so it is shown, not blanked to black.
    inside <- fill_rows_cols(inside)
    brain <- if (is.finite(brain_thr)) {
      brain_component(is.finite(bg) & bg > brain_thr, inside)
    } else inside
    # Orient both images, then lay the checker pattern in display space,
    # anchored at the top-left corner of the visible window, so the key's
    # "top-left tile" statement is true whatever the affine or crop.
    o <- orient_volume_slice_for_raster(bgvol, z, along = along, mat = bg01)
    o_ov <- orient_volume_slice_for_raster(bgvol, z, along = along, mat = ov01)
    o_in <- orient_volume_slice_for_raster(bgvol, z, along = along,
                                           mat = matrix(as.numeric(inside), nrow(bg)))
    o_br <- orient_volume_slice_for_raster(bgvol, z, along = along,
                                           mat = matrix(as.numeric(brain), nrow(bg)))
    x0 <- if (is.null(window)) min(o$x) else window$xlim[1]
    y0 <- if (is.null(window)) max(o$y) else window$ylim[2]
    step <- tile * min(o$display_spacing)
    ci <- floor((o$x - x0) / step)
    ri <- floor((y0 - o$y) / step)
    use_bg <- outer(ri, ci, "+") %% 2 == 0
    chk <- o_ov$mat
    chk[use_bg] <- o$mat[use_bg]
    # Blend from the checkerboard inside the brain to the fixed image alone
    # outside it, over a few voxels, so the mask edge is not a jagged seam.
    # Feather the checker weight over ~3 voxels inside the brain edge, so the
    # boundary is a soft blend rather than a stepped halo.
    wgt <- smooth_mask(o_br$mat > 0, passes = 5L)
    wgt <- pmin(1, pmax(0, (wgt - 0.45) / 0.55))
    chk <- wgt * chk + (1 - wgt) * o$mat
    if (isTRUE(mask_background) && any(o_in$mat > 0)) chk[o_in$mat <= 0] <- NA_real_
    o$mat <- chk
    neuro_tile(o, bg_lim = c(0, 1), bg_cmap = cmap, window = window,
               label = slice_world_label(bgvol, z, along),
               orient = if (isTRUE(annotate) && k == 1L) orientation_letters(o) else NULL,
               tokens = tokens, interpolate = isTRUE(interpolate))
  })
  key <- if (isTRUE(legend)) {
    tile_mm <- tile * min(spacing(space(bgvol))[-along])
    neuro_key_strip(Filter(Negate(is.null), list(
      list(fill = NA, label = sprintf("Top-left tile: %s, alternating with %s", labels[[1L]], labels[[2L]])),
      list(fill = NA, label = sprintf("%s mm tiles; anatomy should run continuously across tile borders",
                                      format_tick(tile_mm))),
      if (isTRUE(focus_brain) && is.finite(brain_thr)) {
        list(fill = NA, label = sprintf("tiles interleave inside the brain; outside it %s only", labels[[1L]]))
      })), tokens)
  }
  o1 <- orient_volume_slice_for_raster(bgvol, zlevels[[1L]], along = along)
  build <- function(cv) neuro_assemble(plots, ncol = ncol, tile_aspect = tile_aspect_of(window, o1),
                        key = key, tokens = tokens, title = title,
                        subtitle = subtitle, caption = caption, canvas = cv)
  cv0 <- canvas_size(canvas)
  fig <- build(cv0)
  if (is.null(canvas)) fig <- neuro_figure(fig, build, cv0, tokens)
  neuro_finish(fig, plots, draw = draw, assemble = assemble, title = title,
               subtitle = subtitle, caption = caption, style = style,
               panel_names = vapply(zlevels, function(z) slice_world_label(bgvol, z, along), ""))
}

#' Box-blur a logical mask into weights in [0, 1]
#' @keywords internal
#' @noRd
smooth_mask <- function(m, passes = 2L) {
  w <- matrix(as.numeric(m), nrow(m), ncol(m))
  w[!is.finite(w)] <- 0
  nr <- nrow(w); nc <- ncol(w)
  if (nr < 3L || nc < 3L) return(w)
  for (k in seq_len(passes)) {
    pad <- matrix(0, nr + 2L, nc + 2L)
    pad[2:(nr + 1L), 2:(nc + 1L)] <- w
    acc <- 0
    for (di in 0:2) for (dj in 0:2) acc <- acc + pad[(1:nr) + di, (1:nc) + dj]
    w <- acc / 9
  }
  w
}
