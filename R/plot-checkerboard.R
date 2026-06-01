#' Checkerboard comparison of two registered volumes
#'
#' Alternates tiles from a background volume and a comparison volume on matched
#' slices. This is useful for visual registration QC.
#'
#' @param bgvol Background/reference 3D volume.
#' @param overlay Comparison 3D volume on the same NeuroSpace grid as `bgvol`.
#' @param zlevels Slices to plot. Defaults to nine evenly spaced slices.
#' @param along Axis for slicing (1 sagittal, 2 coronal, 3 axial).
#' @param tile Tile width in slice pixels.
#' @param cmap Palette used to render the normalized checkerboard image.
#' @param bg_range,ov_range "robust" or "data" intensity scaling.
#' @param probs Quantiles for robust scaling.
#' @param ncol Number of columns in the facet layout.
#' @param title,subtitle,caption Optional layout-level labels used when drawing.
#' @param draw Logical; if `TRUE`, draw the panels on the active graphics
#'   device. If `FALSE`, only return the ggplot objects invisibly.
#' @param style Visual style, either \code{"light"} or \code{"dark"}.
#' @export
plot_checkerboard <- function(
    bgvol, overlay, zlevels = NULL, along = 3L, tile = 16L, cmap = "grays",
    bg_range = c("robust", "data"), ov_range = c("robust", "data"),
    probs = c(.02, .98), ncol = 3L,
    title = NULL, subtitle = NULL, caption = NULL, draw = TRUE,
    style = c("light", "dark")) {
  assert_same_neuro_grid(bgvol, overlay = overlay)
  bg_range <- match.arg(bg_range)
  ov_range <- match.arg(ov_range)
  style <- match.arg(style)
  along <- as.integer(along)
  if (length(along) != 1L || is.na(along) || along < 1L || along > 3L) {
    stop("`along` must be one of 1, 2, or 3.", call. = FALSE)
  }
  tile <- max(1L, as.integer(tile))

  if (is.null(zlevels)) {
    zlevels <- unique(round(seq(1, dim(bgvol)[along], length.out = 9)))
  }
  panel_args <- validate_slice_panel_args(zlevels, along, dim(bgvol), ncol)
  zlevels <- panel_args$zlevels
  along <- panel_args$along
  ncol <- panel_args$ncol

  build_panel <- function(z) {
    sl_bg <- slice(bgvol, z, along = along)
    sl_ov <- slice(overlay, z, along = along)
    bg <- slice_to_matrix(sl_bg)
    ov <- slice_to_matrix(sl_ov)
    if (!identical(dim(bg), dim(ov))) {
      stop("Sliced volumes must have matching dimensions.", call. = FALSE)
    }

    bg_lim <- compute_limits(as.numeric(bg), mode = bg_range, probs = probs)
    ov_lim <- compute_limits(as.numeric(ov), mode = ov_range, probs = probs)
    bg01 <- pmax(0, pmin(1, rescale01(as.numeric(bg), bg_lim)))
    ov01 <- pmax(0, pmin(1, rescale01(as.numeric(ov), ov_lim)))
    bg01 <- matrix(bg01, nrow = nrow(bg), ncol = ncol(bg))
    ov01 <- matrix(ov01, nrow = nrow(ov), ncol = ncol(ov))

    rr <- (row(bg01) - 1L) %/% tile
    cc <- (col(bg01) - 1L) %/% tile
    use_bg <- (rr + cc) %% 2L == 0L
    chk <- ov01
    chk[use_bg] <- bg01[use_bg]

    oriented <- orient_slice_for_raster(sl_bg, chk)
    df <- expand.grid(x = oriented$x, y = oriented$y)
    df$value <- c(t(oriented$mat))
    ggplot2::ggplot(df, ggplot2::aes(x, y, fill = value)) +
      ggplot2::geom_raster(interpolate = FALSE) +
      scale_fill_neuro(cmap = cmap, limits = c(0, 1), guide = "none") +
      ggplot2::coord_fixed() +
      theme_neuro(style = style) +
      ggplot2::labs(title = paste0("z = ", z))
  }

  plots <- lapply(zlevels, build_panel)
  attr(plots, "labels") <- list(title = title, subtitle = subtitle, caption = caption)
  if (!isTRUE(draw)) {
    return(invisible(plots))
  }

  draw_plot_panel_grid(
    plots,
    ncol = ncol,
    title = title,
    subtitle = subtitle,
    caption = caption,
    style = style
  )
}
