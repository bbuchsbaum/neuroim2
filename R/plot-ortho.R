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
  s_ax <- slice(vol, coord[3], along = 3L)  # axial: (i,j)
  s_co <- slice(vol, coord[2], along = 2L)  # coronal: (i,k)
  s_sa <- slice(vol, coord[1], along = 1L)  # sagittal: (j,k)

  # Helper to build world-coordinate data frame for a NeuroSlice
  slice_world_df <- function(sl, downsample) {
    vals <- as.numeric(sl)
    if (downsample > 1L) {
      nr <- dim(sl)[1]; nc <- dim(sl)[2]
      rr <- seq(1L, nr, by = downsample)
      cc <- seq(1L, nc, by = downsample)
      grid <- as.matrix(expand.grid(i = rr, j = cc))
      cds  <- grid_to_coord(space(sl), grid)
      idx  <- grid_to_index(space(sl), grid)
      data.frame(x = cds[,1], y = cds[,2], value = vals[idx])
    } else {
      cds <- index_to_coord(space(sl), seq_len(length(sl)))
      data.frame(x = cds[,1], y = cds[,2], value = vals)
    }
  }

  d_ax <- slice_world_df(s_ax, downsample); d_ax$plane <- "Axial"
  d_co <- slice_world_df(s_co, downsample); d_co$plane <- "Coronal"
  d_sa <- slice_world_df(s_sa, downsample); d_sa$plane <- "Sagittal"

  # Shared limits across panels
  lim <- compute_limits(c(d_ax$value, d_co$value, d_sa$value), range, probs)

  # Crosshair mm positions
  ax_mm <- as.numeric(grid_to_coord(space(s_ax), matrix(c(coord[1], coord[2]), ncol = 2)))
  co_mm <- as.numeric(grid_to_coord(space(s_co), matrix(c(coord[1], coord[3]), ncol = 2)))
  sa_mm <- as.numeric(grid_to_coord(space(s_sa), matrix(c(coord[2], coord[3]), ncol = 2)))

  make_panel_world <- function(df, plane, cross_mm, plane_id) {
    xr <- range(df$x, na.rm = TRUE); yr <- range(df$y, na.rm = TRUE)
    p <- ggplot2::ggplot(df, ggplot2::aes(x, y, fill = value)) +
      ggplot2::geom_raster(interpolate = FALSE) +
      scale_fill_neuro(cmap = cmap, limits = lim) +
      ggplot2::coord_fixed() +
      theme_neuro() +
      ggplot2::labs(title = plane)

    if (crosshair && length(cross_mm) == 2) {
      p <- p +
        ggplot2::geom_segment(x = cross_mm[1], xend = cross_mm[1],
                              y = yr[1], yend = yr[2],
                              linewidth = .3, colour = "white", alpha = .7) +
        ggplot2::geom_segment(y = cross_mm[2], yend = cross_mm[2],
                              x = xr[1], xend = xr[2],
                              linewidth = .3, colour = "white", alpha = .7)
    }
    if (annotate) {
      # Place simple L/R and A/P/S/I based on plane
      labs_tb <- switch(plane_id,
                        axial   = c(top = "A", bottom = "P"),
                        coronal = c(top = "S", bottom = "I"),
                        sagittal= c(top = "S", bottom = "I"))
      p <- p +
        ggplot2::annotate("text", x = xr[1], y = mean(yr), label = "L",
                          colour = "white", fontface = "bold") +
        ggplot2::annotate("text", x = xr[2], y = mean(yr), label = "R",
                          colour = "white", fontface = "bold") +
        ggplot2::annotate("text", x = mean(xr), y = yr[2], label = labs_tb["top"],
                          colour = "white", fontface = "bold") +
        ggplot2::annotate("text", x = mean(xr), y = yr[1], label = labs_tb["bottom"],
                          colour = "white", fontface = "bold")
    }
    p
  }

  pa <- make_panel_world(d_ax, "Axial",    ax_mm, "axial")
  pc <- make_panel_world(d_co, "Coronal",  co_mm, "coronal")
  ps <- make_panel_world(d_sa, "Sagittal", sa_mm, "sagittal")

  # Print in a single row using base grid (no cowplot/patchwork)
  grid::grid.newpage()
  vp <- grid::viewport(layout = grid::grid.layout(1, 3))
  grid::pushViewport(vp)
  print(pa, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(pc, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(ps, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 3))
  invisible(list(axial = pa, coronal = pc, sagittal = ps))
}
