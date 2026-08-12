#' Orthogonal three-plane view with optional crosshairs
#'
#' Creates axial, coronal, and sagittal panels at a given coordinate with
#' harmonized aesthetics. Returns the three ggplot objects invisibly after
#' drawing, or without drawing when \code{draw = FALSE}.
#'
#' @param vol A 3D volume handled by `slice()`.
#' @param coord Length-3 coordinate of the target point. Interpreted as voxel
#'   indices by default; set `unit = "mm"` to convert using `coord_to_grid()`
#'   if available in your environment.
#' @param unit "index" or "mm".
#' @param cmap Palette for the slices.
#' @param range Intensity limits shared by all panels: \code{"robust"},
#'   \code{"data"}, or an explicit numeric \code{c(lo, hi)}.
#' @param probs Quantiles for robust range.
#' @param crosshair Logical; draw crosshair lines.
#' @param annotate Logical; add orientation glyphs.
#' @param downsample Integer decimation for speed.
#' @param title,subtitle,caption Optional layout-level labels used when drawing.
#' @param draw Logical; if `TRUE`, draw the panels on the active graphics
#'   device. If `FALSE`, only return the ggplot objects invisibly.
#' @param style Visual style: \code{"light"}, \code{"dark"}, or \code{"report"}
#'   (light card, dark cropped tiles, typography, and a colorbar -- matching
#'   \code{\link{plot_overlay}}'s report look).
#' @param enhance Display-only enhancement of an unsmoothed statistical
#'   \code{vol}. \code{FALSE} (default) leaves it untouched; \code{TRUE} applies
#'   \code{\link{enhance_stat_map}} with defaults; a named \code{list} is
#'   forwarded as arguments to \code{enhance_stat_map()}.
#' @param crop,interpolate Logical or \code{NULL}; crop panels to the brain
#'   bounding box / smooth the raster. \code{NULL} (default) enables both for
#'   \code{style = "report"} only.
#' @details The affine determines which native voxel axis is nearest each
#'   anatomical plane and how that plane must be permuted or flipped for
#'   display. Oblique images are shown on their regular native voxel planes;
#'   values are not silently resampled. Use \code{deoblique()} or
#'   \code{resample_to()} first when true cardinal-plane sections are required.
#' @export
plot_ortho <- function(
  vol, coord = NULL, unit = c("index","mm"),
  cmap = "grays", range = c("robust","data"), probs = c(.02,.98),
  crosshair = TRUE, annotate = TRUE, downsample = 1L,
  title = NULL, subtitle = NULL, caption = NULL,
  draw = TRUE, style = c("light", "dark", "report"), enhance = FALSE,
  crop = NULL, interpolate = NULL
) {
  unit <- match.arg(unit)
  style <- match.arg(style)

  is_report   <- identical(style, "report")
  panel_style <- .plot_style_colors(style)$panel
  do_crop     <- if (is.null(crop)) is_report else isTRUE(crop)
  interp_bg   <- if (is.null(interpolate)) is_report else isTRUE(interpolate)

  # Optional display-only enhancement of an unsmoothed statistical volume.
  vol <- apply_enhance_arg(vol, enhance)

  if (is.null(coord)) coord <- round(dim(vol) / 2)

  # Convert mm -> voxel grid if possible
  if (unit == "mm") {
    conv <- try(get("coord_to_grid", mode = "function"), silent = TRUE)
    sp   <- try(get("space",        mode = "function"), silent = TRUE)
    if (!inherits(conv, "try-error") && !inherits(sp, "try-error")) {
      coord <- as.integer(conv(sp(vol), matrix(coord, ncol = 3)))
    } else {
      warning("coord_to_grid()/space() not found; treating 'coord' as voxel indices.")
    }
  }
  coord <- as.integer(round(coord))
  if (length(coord) != 3L || anyNA(coord) || any(coord < 1L | coord > dim(vol)[1:3])) {
    stop("`coord` must be a valid length-3 voxel coordinate.", call. = FALSE)
  }

  # Find the native grid axis nearest each anatomical normal. This preserves
  # axial/coronal/sagittal semantics when a valid NIfTI affine permutes the
  # voxel axes, and it keeps oblique native planes on regular raster grids.
  directions <- perm_mat(axes(space(vol)))
  native_axis_for_world <- vapply(
    seq_len(3L),
    function(world_axis) which.max(abs(directions[world_axis, ])),
    integer(1)
  )
  if (anyDuplicated(native_axis_for_world)) {
    stop("Volume axes do not define three distinct anatomical directions.", call. = FALSE)
  }

  make_slice <- function(plane, world_normal) {
    along <- native_axis_for_world[[world_normal]]
    oriented <- orient_volume_slice_for_raster(
      vol,
      z = coord[[along]],
      along = along,
      downsample = downsample
    )
    df <- oriented_raster_df(oriented)
    df$plane <- plane
    list(
      oriented = oriented,
      data = df,
      cross = slice_grid_to_display(oriented, coord[-along])
    )
  }

  axial <- make_slice("Axial", 3L)
  coronal <- make_slice("Coronal", 2L)
  sagittal <- make_slice("Sagittal", 1L)
  d_ax <- axial$data
  d_co <- coronal$data
  d_sa <- sagittal$data

  # Shared limits across panels
  lim <- resolve_display_limits(range, c(d_ax$value, d_co$value, d_sa$value), probs = probs)

  make_panel <- function(slice_info, plane) {
    df <- slice_info$data
    oriented <- slice_info$oriented
    cross_coord <- slice_info$cross
    xr <- range(df$x, na.rm = TRUE); yr <- range(df$y, na.rm = TRUE)

    coord <- ggplot2::coord_fixed()
    if (isTRUE(do_crop)) {
      fin <- df$value[is.finite(df$value)]
      if (length(fin)) {
        thr <- min(fin) + 0.02 * diff(range(fin))
        keep <- is.finite(df$value) & df$value > thr
        if (any(keep)) {
          cx <- range(df$x[keep]); cy <- range(df$y[keep])
          mx <- diff(cx) * 0.06; my <- diff(cy) * 0.06
          coord <- ggplot2::coord_fixed(xlim = c(cx[1] - mx, cx[2] + mx),
                                        ylim = c(cy[1] - my, cy[2] + my),
                                        expand = FALSE)
        }
      }
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x, y, fill = value)) +
      ggplot2::geom_raster(interpolate = interp_bg) +
      scale_fill_neuro(cmap = cmap, limits = lim, guide = "none") +
      coord +
      theme_neuro(style = panel_style) +
      ggplot2::labs(title = plane)
    if (is_report) p <- p + report_tile_theme()

    if (crosshair && length(cross_coord) == 2) {
      p <- p +
        ggplot2::geom_segment(x = cross_coord[1], xend = cross_coord[1],
                              y = yr[1], yend = yr[2],
                              linewidth = .3, colour = "white", alpha = .7) +
        ggplot2::geom_segment(y = cross_coord[2], yend = cross_coord[2],
                              x = xr[1], xend = xr[2],
                              linewidth = .3, colour = "white", alpha = .7)
    }
    if (annotate) {
      p <- p +
        ggplot2::annotate("text", x = xr[1], y = mean(yr), label = oriented$labels[["left"]],
                          colour = "white", fontface = "bold") +
        ggplot2::annotate("text", x = xr[2], y = mean(yr), label = oriented$labels[["right"]],
                          colour = "white", fontface = "bold") +
        ggplot2::annotate("text", x = mean(xr), y = yr[2], label = oriented$labels[["top"]],
                          colour = "white", fontface = "bold") +
        ggplot2::annotate("text", x = mean(xr), y = yr[1], label = oriented$labels[["bottom"]],
                          colour = "white", fontface = "bold")
    }
    p
  }

  pa <- make_panel(axial, "Axial")
  pc <- make_panel(coronal, "Coronal")
  ps <- make_panel(sagittal, "Sagittal")
  plots <- list(axial = pa, coronal = pc, sagittal = ps)
  attr(plots, "labels") <- list(title = title, subtitle = subtitle, caption = caption)

  if (isTRUE(is_report)) {
    combined <- assemble_figure(
      patchwork::wrap_plots(plots, ncol = 3L),
      lim = lim, cmap = cmap, thresh = 0, style = style,
      colorbar = TRUE, title = title, subtitle = subtitle, caption = caption
    )
    if (isTRUE(draw)) print(combined)
    return(invisible(combined))
  }

  if (!isTRUE(draw)) {
    return(invisible(plots))
  }

  draw_plot_panel_grid(
    plots,
    ncol = 3L,
    title = title,
    subtitle = subtitle,
    caption = caption,
    style = style
  )
}
