# Internal utilities used by plotting helpers

#' @keywords internal
#' @noRd
utils::globalVariables(c("x","y","value","z"))

#' Coerce a NeuroSlice (or matrix-like) to a numeric matrix
#' @keywords internal
#' @noRd
slice_to_matrix <- function(slc) {
  if (is.matrix(slc)) {
    return(slc)
  }
  m <- try(as.matrix(slc), silent = TRUE)
  if (!inherits(m, "try-error") && is.matrix(m)) return(m)
  a <- try(as.array(slc), silent = TRUE)
  if (!inherits(a, "try-error")) {
    if (length(dim(a)) == 2) return(a)
    if (length(dim(a)) == 3) return(a[,,1, drop = TRUE])
  }
  # try common slots / accessors without committing to class internals
  data <- try(slc@data, silent = TRUE)
  if (!inherits(data, "try-error") && is.matrix(data)) return(data)
  stop("Cannot coerce 'slc' to a numeric matrix. Provide a 2D matrix or a NeuroSlice with as.matrix().")
}

#' Convert a slice to a tidy data.frame ready for ggplot2::geom_raster
#' @param slc A NeuroSlice or 2D matrix
#' @param downsample Integer decimation factor (>=1)
#' @keywords internal
#' @noRd
slice_df <- function(slc, downsample = 1L) {
  m <- slice_to_matrix(slc)
  if (!is.numeric(m)) m <- suppressWarnings(matrix(as.numeric(m), nrow = nrow(m)))
  if (downsample > 1L) {
    idxr <- seq(1L, nrow(m), by = downsample)
    idxc <- seq(1L, ncol(m), by = downsample)
    m <- m[idxr, idxc, drop = FALSE]
  }
  # Build grid with x=cols, y=rows. y reversed later to match radiological display
  df <- expand.grid(x = seq_len(ncol(m)), y = seq_len(nrow(m)))
  df$value <- c(t(m))
  df
}

#' Compute robust or data-based limits
#' @keywords internal
#' @noRd
compute_limits <- function(x, mode = c("robust","data"), probs = c(.02,.98)) {
  mode <- match.arg(mode)
  x <- x[is.finite(x)]
  if (!length(x)) return(c(0, 1))
  if (mode == "data") {
    rng <- range(x, finite = TRUE)
  } else {
    q <- stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
    rng <- c(min(q), max(q))
    if (!is.finite(rng[1]) || !is.finite(rng[2])) rng <- range(x, finite = TRUE)
    if (rng[1] == rng[2]) rng <- range(x, finite = TRUE)
  }
  rng
}

#' Simple rescale
#' @keywords internal
#' @noRd
rescale01 <- function(x, from) {
  if (is.null(from) || any(!is.finite(from))) return(x)
  a <- from[1]; b <- from[2]
  if (b == a) return(rep(.5, length(x)))
  (x - a) / (b - a)
}

#' Map numeric matrix to RGBA colors
#' @keywords internal
#' @noRd
matrix_to_colors <- function(mat, cmap = "grays", limits = NULL, alpha = 1) {
  cols <- resolve_cmap(cmap, 256)
  if (is.null(limits)) limits <- range(mat[is.finite(mat)])
  s <- rescale01(as.numeric(mat), limits)
  s[!is.finite(s)] <- NA_real_
  # index into palette
  idx <- 1 + floor(s * (length(cols) - 1))
  col <- cols[pmax(1, pmin(length(cols), idx))]
  # apply global alpha
  grDevices::adjustcolor(col, alpha.f = alpha)
}

#' Create a rasterGrob from a numeric matrix using a palette
#' @keywords internal
#' @noRd
matrix_to_raster_grob <- function(mat, cmap = "grays", limits = NULL, alpha = 1) {
  col <- matrix_to_colors(mat, cmap = cmap, limits = limits, alpha = alpha)
  # as.raster expects a matrix of color strings with same dims
  r <- structure(col, dim = dim(mat))
  gr <- grid::rasterGrob(r, interpolate = FALSE)
  gr
}

#' Coordinate helper: fixed aspect and reversed y for radiological convention
#' @keywords internal
#' @noRd
coord_neuro_fixed <- function() {
  list(ggplot2::coord_fixed(), ggplot2::scale_y_reverse(expand = ggplot2::expansion(mult = 0)))
}

#' Add L/R and A/P/S/I annotations (optional)
#'
#' @param plane "axial", "coronal", or "sagittal"
#' @param dims c(nrow, ncol) of the slice matrix
#' @param gp grid::gpar style
#' @return A ggplot2 layer with annotation_custom grobs
#' @export
annotate_orientation <- function(plane = c("axial","coronal","sagittal"),
                                 dims, gp = grid::gpar(col = "white", cex = .9, fontface = "bold")) {
  plane <- match.arg(plane)
  nr <- dims[1]; nc <- dims[2]
  lab_lr <- list(left  = "L", right = "R")
  lab_tb <- switch(plane,
                   axial   = list(top = "A", bottom = "P"),
                   coronal = list(top = "S", bottom = "I"),
                   sagittal= list(top = "S", bottom = "I"))
  layers <- list(
    ggplot2::annotation_custom(grid::textGrob(lab_lr$left,  gp = gp),
                               xmin = 0.5, xmax = 0.5, ymin = nr/2, ymax = nr/2),
    ggplot2::annotation_custom(grid::textGrob(lab_lr$right, gp = gp),
                               xmin = nc + .5, xmax = nc + .5, ymin = nr/2, ymax = nr/2),
    ggplot2::annotation_custom(grid::textGrob(lab_tb$top,   gp = gp),
                               xmin = nc/2, xmax = nc/2, ymin = .5, ymax = .5),
    ggplot2::annotation_custom(grid::textGrob(lab_tb$bottom,gp = gp),
                               xmin = nc/2, xmax = nc/2, ymin = nr + .5, ymax = nr + .5)
  )
  layers
}

