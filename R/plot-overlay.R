#' Composite an overlay map on a structural background
#'
#' Works without extra packages by colorizing both layers to rasters and stacking
#' them as grobs. Great for statistical maps over T1/T2 backgrounds.
#'
#' @param bgvol Background 3D volume.
#' @param overlay Overlay 3D volume on the same NeuroSpace grid as `bgvol`.
#' @param zlevels Slices to plot (indices along the z/3rd axis by default).
#' @param along Axis for slicing (1 sagittal, 2 coronal, 3 axial).
#' @param bg_cmap Background palette (e.g., "grays").
#' @param ov_cmap Overlay palette (e.g., "inferno").
#' @param bg_range,ov_range Background/overlay scaling. Either a mode string,
#'   \code{"robust"} (quantile clip) or \code{"data"} (min/max), or an explicit
#'   numeric \code{c(lo, hi)} to pin the scale (e.g. \code{ov_range = c(-6, 6)})
#'   for consistent coloring across panels and subjects.
#' @param probs Quantiles for robust scaling.
#' @param ov_thresh Numeric threshold; values with |v| < thresh become transparent.
#' @param ov_alpha Global alpha for overlay (0..1).
#' @param ov_alpha_mode One of \code{"binary"} (default: pixels above threshold
#'   get full \code{ov_alpha}, others transparent), \code{"proportional"}
#'   (per-pixel alpha = |v| / cap), or \code{"ramp"} (alpha ramps linearly from 0
#'   at \code{ov_thresh} to 1 at the cap). The cap is shared across all panels so
#'   identical values get identical opacity everywhere.
#' @param ov_symmetric Logical or \code{NULL}. \code{NULL} (default) auto-selects
#'   symmetric limits around zero when the overlay has both positive and negative
#'   values; \code{TRUE}/\code{FALSE} forces the choice. Symmetric limits keep
#'   negative and positive values equally visible with a diverging palette.
#' @param ov_cap Optional numeric; the magnitude used as the upper end of the
#'   (symmetric) color/alpha scale. Defaults to the data-driven limit.
#' @param ncol Number of columns in the panel layout.
#' @param title,subtitle,caption Optional layout-level labels used when drawing.
#' @param draw Logical; if `TRUE`, draw on the active graphics device. If
#'   `FALSE`, return without drawing.
#' @param style Visual style, either \code{"light"} or \code{"dark"}.
#' @param enhance Display-only enhancement of the (unsmoothed) statistical
#'   \code{overlay}. \code{FALSE} (default) leaves it untouched; \code{TRUE}
#'   applies \code{\link{enhance_stat_map}} with defaults; a named \code{list}
#'   is forwarded as arguments to \code{enhance_stat_map()} (e.g.
#'   \code{enhance = list(detail_gain = 2, method = "bilateral")}).
#' @param assemble Logical; if \code{TRUE} (default), return a single assembled
#'   \pkg{patchwork} object (honoring \code{ncol} and the layout labels) suitable
#'   for \code{ggsave()}. If \code{FALSE}, draw a panel grid and return the
#'   per-slice ggplot list invisibly (useful for programmatic access to
#'   individual panels).
#' @param colorbar Logical; when \code{assemble = TRUE}, append a colorbar for
#'   the overlay statistic (with the threshold marked). Default \code{TRUE}.
#' @param legend Logical or \code{NULL}; when \code{assemble = TRUE}, add a
#'   bottom legend strip (positive/negative swatches, threshold, plane). \code{NULL}
#'   (default) shows it for \code{style = "report"} only.
#' @param crop Logical or \code{NULL}; crop every panel to the brain bounding box
#'   (shared across slices, so framing is consistent). \code{NULL} (default)
#'   crops for \code{style = "report"} only.
#' @param interpolate Logical or \code{NULL}; smooth the background raster.
#'   \code{NULL} (default) interpolates for \code{style = "report"} only.
#'
#' @details
#' \strong{Return value.} By default (\code{assemble = TRUE}) the return value is
#' a single \pkg{patchwork} object that can be passed directly to
#' \code{ggsave()}. With \code{assemble = FALSE} the montage is drawn as a side
#' effect and the return value is the \emph{list of per-slice ggplots}
#' (invisibly); passing that list to \code{ggsave()} saves only one panel.
#'
#' \strong{Signed maps.} For overlays with both signs (t/z/contrast maps), the
#' default palette switches to a diverging one and limits become symmetric so
#' negatives are as visible as positives; pass \code{ov_cmap}/\code{ov_symmetric}
#' to override.
#'
#' \strong{Report style.} \code{style = "report"} renders dark brain tiles on a
#' light card with bold/italic typography, a titled colorbar, a bottom legend
#' strip, brain-bbox cropping, and a smoothed background -- a publication-ready
#' look. The individual features (\code{legend}, \code{crop}, \code{interpolate})
#' can also be toggled on any style.
#'
#' @export
plot_overlay <- function(
  bgvol, overlay, zlevels = NULL, along = 3L,
  bg_cmap = "grays", ov_cmap = "inferno",
  bg_range = c("robust","data"), ov_range = c("robust","data"),
  probs = c(.02,.98), ov_thresh = 0, ov_alpha = .7,
  ov_alpha_mode = c("binary", "proportional", "ramp"), ov_symmetric = NULL,
  ov_cap = NULL, ncol = 3L, title = NULL, subtitle = NULL, caption = NULL,
  draw = TRUE, style = c("light", "dark", "report"), enhance = FALSE,
  assemble = TRUE, colorbar = TRUE, legend = NULL, crop = NULL, interpolate = NULL
) {
  ov_cmap_missing <- missing(ov_cmap)
  assert_same_neuro_grid(bgvol, overlay = overlay)
  ov_alpha_mode <- match.arg(ov_alpha_mode)
  style <- match.arg(style)

  # Per-style feature defaults (each independently overridable).
  is_report   <- identical(style, "report")
  panel_style <- .plot_style_colors(style)$panel
  show_legend <- if (is.null(legend)) is_report else isTRUE(legend)
  do_crop     <- if (is.null(crop)) is_report else isTRUE(crop)
  interp_bg   <- if (is.null(interpolate)) is_report else isTRUE(interpolate)

  # Optional display-only enhancement of the (unsmoothed) statistical overlay.
  overlay <- apply_enhance_arg(overlay, enhance)
  along <- as.integer(along)
  if (length(along) != 1L || is.na(along) || along < 1L || along > 3L) {
    stop("`along` must be one of 1, 2, or 3.", call. = FALSE)
  }

  if (is.null(zlevels)) zlevels <- unique(round(seq(1, dim(bgvol)[along], length.out = 9)))
  panel_args <- validate_slice_panel_args(zlevels, along, dim(bgvol), ncol)
  zlevels <- panel_args$zlevels
  along <- panel_args$along
  ncol <- panel_args$ncol

  selected_values <- function(vol) {
    unlist(lapply(zlevels, function(z) as.numeric(slice_to_matrix(slice(vol, z, along = along)))))
  }
  ov_vals <- selected_values(overlay)
  bg_lim <- resolve_display_limits(bg_range, selected_values(bgvol), probs = probs)
  ov_lim <- resolve_display_limits(ov_range, ov_vals, probs = probs)

  # Detect a signed (diverging) statistical map across the displayed slices.
  finite_ov <- ov_vals[is.finite(ov_vals)]
  signed <- length(finite_ov) > 0L && min(finite_ov) < 0 && max(finite_ov) > 0
  symmetric <- if (is.null(ov_symmetric)) signed else isTRUE(ov_symmetric)

  # Choose an appropriate palette for signed data. When the caller did not pick a
  # palette we silently default to a diverging one (documented); when they
  # explicitly chose a sequential palette for signed data we warn, since
  # negatives will render on the dark end and read as holes.
  if (signed) {
    if (ov_cmap_missing) {
      ov_cmap <- "blue-red"
    } else if (!is_diverging_cmap(ov_cmap)) {
      warning("Overlay has both positive and negative values but `ov_cmap` is not a ",
              "diverging palette; negative values may render poorly. Consider a ",
              "diverging palette such as 'RdBu' or 'blue-red'.", call. = FALSE)
    }
  }

  # Shared color/alpha cap so identical values look identical across all panels.
  cap <- if (!is.null(ov_cap)) abs(ov_cap) else max(abs(ov_lim), na.rm = TRUE)
  if (symmetric && is.finite(cap) && cap > 0) {
    ov_lim <- c(-cap, cap)
  }
  ov_abs_max <- max(abs(ov_lim), na.rm = TRUE)
  if (!is.finite(ov_abs_max) || ov_abs_max <= 0) ov_abs_max <- 1

  # Shared brain bounding box (so every panel is framed identically).
  crop_win <- NULL
  if (isTRUE(do_crop)) {
    bg_fin <- selected_values(bgvol)
    bg_fin <- bg_fin[is.finite(bg_fin)]
    if (length(bg_fin)) {
      bg_thresh <- min(bg_fin) + 0.02 * diff(range(bg_fin))
      crop_win <- compute_crop_window(bgvol, overlay, zlevels, along,
                                      bg_thresh = bg_thresh, ov_thresh = ov_thresh)
    }
  }

  build_panel <- function(z) {
    sl_bg <- slice(bgvol, z, along = along)
    sl_ov <- slice(overlay, z, along = along)

    bg <- slice_to_matrix(sl_bg)
    bg_oriented <- orient_slice_for_raster(sl_bg, bg)
    df_bg <- expand.grid(x = bg_oriented$x, y = bg_oriented$y)
    df_bg$value <- c(t(bg_oriented$mat))

    # Overlay data; apply threshold then convert to grob for independent palette
    mov <- slice_to_matrix(sl_ov)
    if (ov_alpha_mode == "binary") {
      if (isTRUE(ov_thresh > 0)) {
        below <- !is.na(mov) & abs(mov) < ov_thresh
        mov[below] <- NA_real_
      }
      oriented <- orient_slice_for_raster(sl_ov, mov)
      g_ov <- matrix_to_raster_grob(
        oriented$mat,
        cmap = ov_cmap,
        limits = ov_lim,
        alpha = ov_alpha
      )
    } else {
      abs_mov <- abs(mov)
      if (ov_alpha_mode == "proportional") {
        amap <- abs_mov / ov_abs_max            # shared (not per-slice) denominator
      } else {                                   # "ramp": 0 at thresh -> 1 at cap
        denom <- ov_abs_max - ov_thresh
        if (!is.finite(denom) || denom <= 0) denom <- ov_abs_max
        amap <- (abs_mov - ov_thresh) / denom
      }
      amap[] <- pmin(pmax(amap, 0), 1)
      if (isTRUE(ov_thresh > 0)) {
        amap[!is.na(abs_mov) & abs_mov < ov_thresh] <- 0
      }
      amap[is.na(abs_mov)] <- 0
      oriented <- orient_slice_for_raster(sl_ov, mov, alpha_map = amap)
      g_ov <- matrix_to_raster_grob(
        oriented$mat,
        cmap = ov_cmap,
        limits = ov_lim,
        alpha = ov_alpha,
        alpha_map = oriented$alpha_map
      )
    }

    slice_label <- switch(as.character(along), "1" = "x", "2" = "y", "3" = "z")
    xr <- raster_extent_from_centers(bg_oriented$x)
    yr <- raster_extent_from_centers(bg_oriented$y)

    coord <- if (is.null(crop_win)) {
      ggplot2::coord_fixed()
    } else {
      ggplot2::coord_fixed(xlim = crop_win$xlim, ylim = crop_win$ylim, expand = FALSE)
    }

    p <- ggplot2::ggplot(df_bg, ggplot2::aes(x, y, fill = value)) +
      ggplot2::geom_raster(interpolate = interp_bg) +
      scale_fill_neuro(cmap = bg_cmap, limits = bg_lim, guide = "none") +
      coord +
      theme_neuro(style = panel_style) +
      ggplot2::labs(title = paste0(slice_label, " = ", z)) +
      ggplot2::annotation_custom(g_ov, xmin = xr[1], xmax = xr[2], ymin = yr[1], ymax = yr[2])

    if (is_report) {
      # Borderless dark tiles with a touch more gutter, to read as cards on the
      # light background.
      p <- p + ggplot2::theme(
        panel.border = ggplot2::element_blank(),
        plot.margin  = grid::unit(c(5, 5, 5, 5), "pt")
      )
    }

    p
  }

  plots <- lapply(zlevels, build_panel)
  attr(plots, "labels") <- list(title = title, subtitle = subtitle, caption = caption)

  if (isTRUE(assemble)) {
    top <- patchwork::wrap_plots(plots, ncol = ncol)
    if (isTRUE(colorbar)) {
      cbar <- make_overlay_colorbar(ov_lim, ov_cmap, thresh = ov_thresh,
                                    style = style, title = "value")
      top <- patchwork::wrap_plots(top, cbar, nrow = 1L, widths = c(1, 0.08))
    }

    combined <- top
    if (isTRUE(show_legend)) {
      pal <- resolve_cmap(ov_cmap, 256)
      plane <- switch(as.character(along), "1" = "Sagittal", "2" = "Coronal", "Axial")
      leg <- make_overlay_legend(ov_thresh, pal[length(pal)], pal[1L],
                                 symmetric = symmetric, style = style, plane = plane)
      combined <- patchwork::wrap_plots(top, leg, ncol = 1L, heights = c(1, 0.11))
    }

    combined <- combined +
      patchwork::plot_annotation(title = title, subtitle = subtitle, caption = caption,
                                 theme = annotation_theme(style))
    if (isTRUE(draw)) print(combined)
    return(invisible(combined))
  }

  if (!isTRUE(draw)) {
    return(invisible(plots))
  }

  draw_plot_panel_grid(
    plots,
    ncol = ncol,
    title = title,
    subtitle = subtitle,
    caption = caption,
    style = if (is_report) "dark" else style
  )
}
