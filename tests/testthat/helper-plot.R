# --- inspection helpers for the plot_* family --------------------------------
#
# A slice tile is a ggplot whose layers are: the background geom_raster, then
# one annotation_custom per overlay raster, then text annotations (slice label,
# orientation letters), then (plot_ortho only) crosshair segments. These
# helpers locate layers by what they contain rather than by fixed position.

#' Collect every text label drawn by a (possibly assembled) figure.
grob_labels <- function(p) {
  collect <- function(g) {
    c(if (inherits(g, "text")) as.character(g$label),
      unlist(lapply(g$grobs, collect)),
      unlist(lapply(g$children, collect)))
  }
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  unname(collect(grid::grid.grabExpr(print(p))))
}

layer_geom <- function(layer) class(layer$geom)[[1L]]

#' Built data of every layer whose geom is `geom` (e.g. "GeomText").
built_layers_of <- function(p, geom) {
  built <- ggplot2::ggplot_build(p)
  keep <- vapply(p$layers, layer_geom, character(1)) == geom
  built$data[keep]
}

#' Built data of the foreground text layers of a tile.
#'
#' Tile text is drawn over a translucent "halo" (offset copies of the same
#' labels); the halo layers carry alpha < 1 and are skipped here.
tile_text_data <- function(p) {
  Filter(function(d) !("alpha" %in% names(d)) || all(is.na(d$alpha) | d$alpha >= 1),
         built_layers_of(p, "GeomText"))
}

#' All text annotations of a tile, in drawing order.
tile_text_labels <- function(p) {
  unlist(lapply(tile_text_data(p), function(d) as.character(d$label)))
}

#' The in-tile slice label (the first text annotation).
tile_slice_label <- function(p) tile_text_labels(p)[[1L]]

#' Orientation letters of a tile, named by the side of the tile they sit on.
#'
#' Letters are single-character text annotations; their side is read from
#' their position relative to the centre of the tile's coordinate window.
tile_orientation_letters <- function(p) {
  lim <- p$coordinates$limits
  cx <- mean(lim$x)
  cy <- mean(lim$y)
  out <- character()
  for (d in tile_text_data(p)) {
    for (i in seq_len(nrow(d))) {
      lab <- as.character(d$label[[i]])
      if (nchar(lab) != 1L) next
      dx <- d$x[[i]] - cx
      dy <- d$y[[i]] - cy
      side <- if (abs(dx) > abs(dy)) {
        if (dx < 0) "left" else "right"
      } else {
        if (dy > 0) "top" else "bottom"
      }
      out[[side]] <- lab
    }
  }
  out
}

#' Crosshair segments of an ortho tile as a data frame (x, xend, y, yend).
tile_crosshair <- function(p) {
  segs <- built_layers_of(p, "GeomSegment")
  if (!length(segs)) return(NULL)
  do.call(rbind, lapply(segs, function(d) d[, c("x", "xend", "y", "yend")]))
}

#' A world-coordinate slice label such as "z = -12" (optionally suffixed " mm").
expect_world_label <- function(label, expected) {
  expect_match(label, paste0("^", expected, "( mm)?$"))
}
