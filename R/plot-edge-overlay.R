#' Overlay fixed and moving edge maps on a background volume
#'
#' Displays a structural/reference background with two edge channels rendered in
#' distinct colors. This is intended for registration QC where fixed/template
#' edges and moving/result edges should coincide.
#'
#' @param bgvol Background 3D volume.
#' @param fixed_edges Edge map for the fixed/reference image on the same
#'   NeuroSpace grid as `bgvol`.
#' @param moving_edges Edge map for the moving/result image on the same
#'   NeuroSpace grid as `bgvol`.
#' @param zlevels Slices to plot. Defaults to nine evenly spaced slices.
#' @param along Axis for slicing (1 sagittal, 2 coronal, 3 axial).
#' @param bg_cmap Background palette.
#' @param fixed_color,moving_color Overlay colors for the two edge maps.
#' @param bg_range,edge_range "robust" or "data" intensity scaling.
#' @param probs Quantiles for robust scaling.
#' @param edge_thresh Values below this edge magnitude are transparent.
#' @param edge_alpha Global alpha for edge overlays.
#' @param ncol Number of columns in the facet layout.
#' @param title,subtitle,caption Optional layout-level labels used when drawing.
#' @param draw Logical; if `TRUE`, draw the panels on the active graphics
#'   device. If `FALSE`, only return the ggplot objects invisibly.
#' @param style Visual style, either \code{"light"} or \code{"dark"}.
#' @export
plot_edge_overlay <- function(
    bgvol, fixed_edges, moving_edges, zlevels = NULL, along = 3L,
    bg_cmap = "grays", fixed_color = "#00d5ff", moving_color = "#ff3b30",
    bg_range = c("robust", "data"), edge_range = c("robust", "data"),
    probs = c(.02, .98), edge_thresh = 0, edge_alpha = .85,
    ncol = 3L, title = NULL, subtitle = NULL, caption = NULL, draw = TRUE,
    style = c("light", "dark")) {
  assert_same_neuro_grid(bgvol, fixed_edges = fixed_edges, moving_edges = moving_edges)
  bg_range <- match.arg(bg_range)
  edge_range <- match.arg(edge_range)
  style <- match.arg(style)
  along <- as.integer(along)
  if (length(along) != 1L || is.na(along) || along < 1L || along > 3L) {
    stop("`along` must be one of 1, 2, or 3.", call. = FALSE)
  }

  if (is.null(zlevels)) {
    zlevels <- unique(round(seq(1, dim(bgvol)[along], length.out = 9)))
  }
  panel_args <- validate_slice_panel_args(zlevels, along, dim(bgvol), ncol)
  zlevels <- panel_args$zlevels
  along <- panel_args$along
  ncol <- panel_args$ncol

  edge_alpha_map <- function(mat, limits) {
    vals <- abs(mat)
    if (any(!is.finite(limits)) || limits[[1L]] == limits[[2L]]) {
      alpha <- as.numeric(vals > edge_thresh & vals > 0)
    } else {
      alpha <- pmax(0, pmin(1, rescale01(as.numeric(vals), limits)))
    }
    alpha <- matrix(alpha, nrow = nrow(mat), ncol = ncol(mat))
    alpha[!is.finite(alpha)] <- 0
    alpha[vals <= edge_thresh] <- 0
    alpha
  }

  build_panel <- function(z) {
    sl_bg <- slice(bgvol, z, along = along)
    sl_fix <- slice(fixed_edges, z, along = along)
    sl_mov <- slice(moving_edges, z, along = along)

    bg <- slice_to_matrix(sl_bg)
    fix <- abs(slice_to_matrix(sl_fix))
    mov <- abs(slice_to_matrix(sl_mov))
    bg_oriented <- orient_slice_for_raster(sl_bg, bg)
    df <- expand.grid(x = bg_oriented$x, y = bg_oriented$y)
    df$value <- c(t(bg_oriented$mat))

    bg_lim <- compute_limits(as.numeric(bg), mode = bg_range, probs = probs)
    edge_lim <- compute_limits(c(as.numeric(fix), as.numeric(mov)), mode = edge_range, probs = probs)

    fix_alpha <- edge_alpha_map(fix, edge_lim)
    mov_alpha <- edge_alpha_map(mov, edge_lim)
    fix_oriented <- orient_slice_for_raster(sl_fix, fix, alpha_map = fix_alpha)
    mov_oriented <- orient_slice_for_raster(sl_mov, mov, alpha_map = mov_alpha)

    g_fix <- matrix_to_raster_grob(
      fix_oriented$mat,
      cmap = c(fixed_color, fixed_color),
      limits = edge_lim,
      alpha = edge_alpha,
      alpha_map = fix_oriented$alpha_map
    )
    g_mov <- matrix_to_raster_grob(
      mov_oriented$mat,
      cmap = c(moving_color, moving_color),
      limits = edge_lim,
      alpha = edge_alpha,
      alpha_map = mov_oriented$alpha_map
    )

    xr <- raster_extent_from_centers(bg_oriented$x)
    yr <- raster_extent_from_centers(bg_oriented$y)
    ggplot2::ggplot(df, ggplot2::aes(x, y, fill = value)) +
      ggplot2::geom_raster(interpolate = FALSE) +
      scale_fill_neuro(cmap = bg_cmap, limits = bg_lim, guide = "none") +
      ggplot2::coord_fixed() +
      theme_neuro(style = style) +
      ggplot2::labs(title = paste0("z = ", z)) +
      ggplot2::annotation_custom(g_fix, xmin = xr[1], xmax = xr[2], ymin = yr[1], ymax = yr[2]) +
      ggplot2::annotation_custom(g_mov, xmin = xr[1], xmax = xr[2], ymin = yr[1], ymax = yr[2])
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
