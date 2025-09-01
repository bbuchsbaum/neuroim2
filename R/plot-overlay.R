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
    mbg <- slice_to_matrix(sl_bg)
    mov <- slice_to_matrix(sl_ov)

    # Compute limits per layer
    bg_lim <- compute_limits(as.numeric(mbg), mode = bg_range, probs = probs)
    ov_lim <- compute_limits(as.numeric(mov), mode = ov_range, probs = probs)

    # Threshold overlay
    if (isTRUE(ov_thresh > 0)) {
      mov[abs(mov) < ov_thresh] <- NA_real_
    }

    # Build grobs
    g_bg <- matrix_to_raster_grob(mbg, cmap = bg_cmap, limits = bg_lim, alpha = 1)
    g_ov <- matrix_to_raster_grob(mov, cmap = ov_cmap, limits = ov_lim, alpha = ov_alpha)

    nr <- nrow(mbg); nc <- ncol(mbg)
    df_stub <- data.frame(x = c(1, nc), y = c(1, nr), value = c(0, 1))

    # Compose using annotation_custom with explicit extents
    p <- ggplot2::ggplot(df_stub, ggplot2::aes(x, y)) +
      ggplot2::geom_blank() +
      ggplot2::annotation_custom(g_bg, xmin = .5, xmax = nc + .5, ymin = .5, ymax = nr + .5) +
      ggplot2::annotation_custom(g_ov, xmin = .5, xmax = nc + .5, ymin = .5, ymax = nr + .5) +
      coord_neuro_fixed() +
      theme_neuro() +
      ggplot2::labs(title = paste0("z = ", z))

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

