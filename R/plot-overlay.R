#' Composite an overlay map on a structural background
#'
#' Works without extra packages by colorizing both layers to rasters and stacking
#' them as grobs. Great for statistical maps over T1/T2 backgrounds.
#'
#' @param bgvol Background 3D volume.
#' @param overlay Overlay 3D volume (same dims as bgvol).
#' @param zlevels Slices to plot (indices along the z/3rd axis by default).
#' @param along Axis for slicing (1 sagittal, 2 coronal, 3 axial).
#' @param bg_cmap Background palette (e.g., "grays").
#' @param ov_cmap Overlay palette (e.g., "inferno").
#' @param bg_range,ov_range "robust" or "data" for background/overlay scaling.
#' @param probs Quantiles for robust scaling.
#' @param ov_thresh Numeric threshold; values with |v| < thresh become transparent.
#' @param ov_alpha Global alpha for overlay (0..1).
#' @param ncol Number of columns in the facet layout.
#' @param title,subtitle,caption Optional labels.
#' @export
plot_overlay <- function(
  bgvol, overlay, zlevels = NULL, along = 3L,
  bg_cmap = "grays", ov_cmap = "inferno",
  bg_range = c("robust","data"), ov_range = c("robust","data"),
  probs = c(.02,.98), ov_thresh = 0, ov_alpha = .7,
  ncol = 3L, title = NULL, subtitle = NULL, caption = NULL
) {
  stopifnot(all(dim(bgvol) == dim(overlay)))
  bg_range <- match.arg(bg_range); ov_range <- match.arg(ov_range)

  if (is.null(zlevels)) zlevels <- unique(round(seq(1, dim(bgvol)[along], length.out = 6)))

  build_panel <- function(z) {
    sl_bg <- slice(bgvol, z, along = along)
    sl_ov <- slice(overlay, z, along = along)

    # Background data in world (mm) coordinates
    vals_bg <- as.numeric(sl_bg)
    cds_bg  <- index_to_coord(space(sl_bg), seq_len(length(sl_bg)))
    df_bg   <- data.frame(x = cds_bg[,1], y = cds_bg[,2], value = vals_bg)

    # Compute limits per layer
    bg_lim <- compute_limits(df_bg$value, mode = bg_range, probs = probs)

    # Overlay data; apply threshold then convert to grob for independent palette
    mov <- slice_to_matrix(sl_ov)
    if (isTRUE(ov_thresh > 0)) mov[abs(mov) < ov_thresh] <- NA_real_
    ov_lim <- compute_limits(as.numeric(mov), mode = ov_range, probs = probs)
    g_ov <- matrix_to_raster_grob(mov, cmap = ov_cmap, limits = ov_lim, alpha = ov_alpha)

    xr <- range(df_bg$x, na.rm = TRUE)
    yr <- range(df_bg$y, na.rm = TRUE)

    p <- ggplot2::ggplot(df_bg, ggplot2::aes(x, y, fill = value)) +
      ggplot2::geom_raster(interpolate = FALSE) +
      scale_fill_neuro(cmap = bg_cmap, limits = bg_lim) +
      ggplot2::coord_fixed() +
      theme_neuro() +
      ggplot2::labs(title = paste0("z = ", z)) +
      ggplot2::annotation_custom(g_ov, xmin = xr[1], xmax = xr[2], ymin = yr[1], ymax = yr[2])

    p
  }

  plots <- lapply(zlevels, build_panel)

  # Arrange via grid viewports
  grid::grid.newpage()
  n <- length(plots)
  ncol <- min(ncol, n)
  nrow <- ceiling(n / ncol)
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow, ncol)))
  for (i in seq_len(n)) {
    r <- ((i - 1) %/% ncol) + 1
    c <- ((i - 1) %%  ncol) + 1
    print(plots[[i]], vp = grid::viewport(layout.pos.row = r, layout.pos.col = c))
  }
  grid::upViewport(0)
  invisible(plots)
}
