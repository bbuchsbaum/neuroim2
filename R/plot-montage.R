#' Plot a montage of axial (or any-plane) slices using facetting
#'
#' This avoids extra dependencies by using a single ggplot with facets
#' and a shared colorbar. Supply a list of slice objects or a volume + indices.
#'
#' @param x Either a 3D volume object accepted by `slice()` or a list of slices.
#' @param zlevels Integer indices of slices to plot (if `x` is a volume).
#' @param along Axis along which to slice (1 = sagittal, 2 = coronal, 3 = axial).
#' @param cmap Palette name or vector (see [resolve_cmap()]).
#' @param range "robust" (quantile-based), "data" (min/max), or an explicit
#'   numeric \code{c(lo, hi)}.
#' @param probs Quantiles for `range="robust"`.
#' @param ncol Number of columns in the facet layout.
#' @param downsample Integer decimation for speed.
#' @param title,subtitle,caption Optional ggplot labels.
#' @param style Visual style: \code{"light"}, \code{"dark"}, or \code{"report"}
#'   (light card, dark cropped tiles, typography, and a colorbar -- matching
#'   \code{\link{plot_overlay}}'s report look).
#' @param crop,interpolate Logical or \code{NULL}; crop to the brain bounding box
#'   / smooth the raster. \code{NULL} (default) enables both for
#'   \code{style = "report"} only. (Cropping applies to the volume path.)
#' @export
plot_montage <- function(
  x, zlevels = NULL, along = 3L,
  cmap = "grays", range = c("robust","data"), probs = c(.02,.98),
  ncol = 6L, downsample = 1L,
  title = NULL, subtitle = NULL, caption = NULL,
  style = c("light", "dark", "report"),
  crop = NULL, interpolate = NULL
) {
  style <- match.arg(style)
  is_report <- identical(style, "report")
  do_crop   <- if (is.null(crop)) is_report else isTRUE(crop)
  interp_bg <- if (is.null(interpolate)) is_report else isTRUE(interpolate)

  world <- inherits(x, "NeuroVol")
  if (world) {
    # World-coordinate aware path (respects anatomical orientation)
    if (is.null(zlevels)) {
      zlevels <- unique(round(seq(1, dim(x)[along], length.out = 12)))
    }
    dfl <- lapply(zlevels, function(z) {
      sl <- slice(x, z, along = along)
      vals <- as.numeric(sl)
      # Optional downsample by grid decimation
      if (downsample > 1L) {
        nr <- dim(sl)[1]; nc <- dim(sl)[2]
        rr <- seq(1L, nr, by = downsample)
        cc <- seq(1L, nc, by = downsample)
        grid <- as.matrix(expand.grid(i = rr, j = cc))
        cds  <- grid_to_coord(space(sl), grid)
        idx  <- grid_to_index(space(sl), grid)
        data.frame(x = cds[,1], y = cds[,2], z = z, value = vals[idx])
      } else {
        cds <- index_to_coord(space(sl), seq_len(length(sl)))
        data.frame(x = cds[,1], y = cds[,2], z = z, value = vals)
      }
    })
  } else {
    # Fallback for lists of slices / plain matrices (pixel grid)
    make_slice <- function(z) {
      sl <- slice(x, z, along = along)
      df <- slice_df(sl, downsample = downsample)
      df$z <- z
      df
    }
    if (is.list(x) && is.null(zlevels)) {
      dfl <- Map(function(sl, idx) { df <- slice_df(sl, downsample); df$z <- idx; df },
                 x, seq_along(x))
    } else {
      if (is.null(zlevels)) {
        zlevels <- unique(round(seq(1, dim(x)[along], length.out = 12)))
      }
      dfl <- lapply(zlevels, make_slice)
    }
  }
  df <- do.call(rbind, dfl)
  lim <- resolve_display_limits(range, df$value, probs = probs)

  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, fill = value)) +
    ggplot2::geom_raster(interpolate = interp_bg) +
    ggplot2::facet_wrap(~ z, ncol = ncol, scales = "fixed") +
    scale_fill_neuro(cmap = cmap, limits = lim,
                     guide = if (is_report) "none" else "colourbar")

  # Crop to the brain bounding box (volume/world path only).
  crop_coord <- NULL
  if (isTRUE(do_crop) && world) {
    fin <- df$value[is.finite(df$value)]
    if (length(fin)) {
      thr  <- min(fin) + 0.02 * diff(range(fin))
      keep <- is.finite(df$value) & df$value > thr
      if (any(keep)) {
        cx <- range(df$x[keep]); cy <- range(df$y[keep])
        mx <- diff(cx) * 0.06;   my <- diff(cy) * 0.06
        crop_coord <- ggplot2::coord_fixed(xlim = c(cx[1] - mx, cx[2] + mx),
                                           ylim = c(cy[1] - my, cy[2] + my),
                                           expand = FALSE)
      }
    }
  }
  p <- p + if (!is.null(crop_coord)) crop_coord
           else if (world) ggplot2::coord_fixed() else coord_neuro_fixed()

  if (is_report) {
    p <- p + report_facet_theme()
    return(assemble_figure(p, lim = lim, cmap = cmap, thresh = 0, style = "report",
                           colorbar = TRUE, title = title, subtitle = subtitle,
                           caption = caption))
  }

  p +
    theme_neuro(style = style) +
    ggplot2::labs(title = title, subtitle = subtitle, caption = caption)
}
