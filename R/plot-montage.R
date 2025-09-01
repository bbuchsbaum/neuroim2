#' Plot a montage of axial (or any-plane) slices using facetting
#'
#' This avoids extra dependencies by using a single ggplot with facets
#' and a shared colorbar. Supply a list of slice objects or a volume + indices.
#'
#' @param x Either a 3D volume object accepted by `slice()` or a list of slices.
#' @param zlevels Integer indices of slices to plot (if `x` is a volume).
#' @param along Axis along which to slice (1 = sagittal, 2 = coronal, 3 = axial).
#' @param cmap Palette name or vector (see [resolve_cmap()]).
#' @param range "robust" (quantile-based) or "data" (min/max).
#' @param probs Quantiles for `range="robust"`.
#' @param ncol Number of columns in the facet layout.
#' @param downsample Integer decimation for speed.
#' @param title,subtitle,caption Optional ggplot labels.
#' @export
plot_montage <- function(
  x, zlevels = NULL, along = 3L,
  cmap = "grays", range = c("robust","data"), probs = c(.02,.98),
  ncol = 6L, downsample = 1L,
  title = NULL, subtitle = NULL, caption = NULL
) {
  range <- match.arg(range)
  if (inherits(x, "NeuroVol")) {
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
    df <- do.call(rbind, dfl)
    lim <- compute_limits(df$value, mode = range, probs = probs)

    p <- ggplot2::ggplot(df, ggplot2::aes(x, y, fill = value)) +
      ggplot2::geom_raster(interpolate = FALSE) +
      ggplot2::facet_wrap(~ z, ncol = ncol, scales = "fixed") +
      scale_fill_neuro(cmap = cmap, limits = lim) +
      ggplot2::coord_fixed() +
      theme_neuro() +
      ggplot2::labs(title = title, subtitle = subtitle, caption = caption)
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
    df <- do.call(rbind, dfl)
    lim <- compute_limits(df$value, mode = range, probs = probs)

    p <- ggplot2::ggplot(df, ggplot2::aes(x, y, fill = value)) +
      ggplot2::geom_raster(interpolate = FALSE) +
      ggplot2::facet_wrap(~ z, ncol = ncol, scales = "fixed") +
      scale_fill_neuro(cmap = cmap, limits = lim) +
      coord_neuro_fixed() +
      theme_neuro() +
      ggplot2::labs(title = title, subtitle = subtitle, caption = caption)
  }

  p
}
