#' Orthogonal three-plane view with optional crosshairs
#'
#' Creates axial, coronal, and sagittal panels at a given coordinate with
#' harmonized aesthetics. Returns (invisibly) the three ggplot objects after
#' printing them in a single row using base grid (no extra deps).
#'
#' @param vol A 3D volume handled by `slice()`.
#' @param coord Length-3 coordinate of the target point. Interpreted as voxel
#'   indices by default; set `unit = "mm"` to convert using `coord_to_grid()`
#'   if available in your environment.
#' @param unit "index" or "mm".
#' @param cmap Palette for the slices.
#' @param range "robust" or "data" for intensity limits shared by all panels.
#' @param probs Quantiles for robust range.
#' @param crosshair Logical; draw crosshair lines.
#' @param annotate Logical; add orientation glyphs.
#' @param downsample Integer decimation for speed.
#' @export
plot_ortho <- function(
  vol, coord = NULL, unit = c("index","mm"),
  cmap = "grays", range = c("robust","data"), probs = c(.02,.98),
  crosshair = TRUE, annotate = TRUE, downsample = 1L
) {
  unit <- match.arg(unit)
  range <- match.arg(range)

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

  # Extract three orthogonal slices
  s_ax <- slice(vol, coord[3], along = 3L)
  s_co <- slice(vol, coord[2], along = 2L)
  s_sa <- slice(vol, coord[1], along = 1L)

  # Build data frames
  d_ax <- slice_df(s_ax, downsample); d_ax$plane <- "Axial"
  d_co <- slice_df(s_co, downsample); d_co$plane <- "Coronal"
  d_sa <- slice_df(s_sa, downsample); d_sa$plane <- "Sagittal"

  # Shared limits
  lim <- compute_limits(c(d_ax$value, d_co$value, d_sa$value), range, probs)

  make_panel <- function(df, plane, crosshair_xy = NULL, dims, plane_id) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x, y, fill = value)) +
      ggplot2::geom_raster(interpolate = FALSE) +
      scale_fill_neuro(cmap = cmap, limits = lim) +
      coord_neuro_fixed() +
      theme_neuro() +
      ggplot2::labs(title = plane)

    if (crosshair && !is.null(crosshair_xy)) {
      p <- p +
        ggplot2::geom_segment(x = crosshair_xy$x, xend = crosshair_xy$x,
                              y = 0.5, yend = dims[1] + .5,
                              linewidth = .3, colour = "white", alpha = .7) +
        ggplot2::geom_segment(y = crosshair_xy$y, yend = crosshair_xy$y,
                              x = 0.5, xend = dims[2] + .5,
                              linewidth = .3, colour = "white", alpha = .7)
    }
    if (annotate) {
      p <- p + annotate_orientation(plane = plane_id, dims = dims)
    }
    p
  }

  # For crosshair positions, we use the voxel indices in the respective plane
  # Axial (z = coord[3]): x->coord[1], y->coord[2]
  ma <- slice_to_matrix(s_ax); pa <- make_panel(d_ax, "Axial",
         list(x = coord[1] / downsample, y = coord[2] / downsample),
         dims = c(nrow(ma), ncol(ma)), plane_id = "axial")
  # Coronal (y = coord[2]): x->coord[1], y->coord[3]
  mc <- slice_to_matrix(s_co); pc <- make_panel(d_co, "Coronal",
         list(x = coord[1] / downsample, y = coord[3] / downsample),
         dims = c(nrow(mc), ncol(mc)), plane_id = "coronal")
  # Sagittal (x = coord[1]): x->coord[2], y->coord[3]
  ms <- slice_to_matrix(s_sa); ps <- make_panel(d_sa, "Sagittal",
         list(x = coord[2] / downsample, y = coord[3] / downsample),
         dims = c(nrow(ms), ncol(ms)), plane_id = "sagittal")

  # Print in a single row using base grid (no cowplot/patchwork)
  grid::grid.newpage()
  vp <- grid::viewport(layout = grid::grid.layout(1, 3))
  grid::pushViewport(vp)
  print(pa, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(pc, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(ps, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 3))
  invisible(list(axial = pa, coronal = pc, sagittal = ps))
}

