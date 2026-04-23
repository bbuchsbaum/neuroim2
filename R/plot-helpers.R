# Internal utilities used by plotting helpers

#' @keywords internal
#' @noRd
utils::globalVariables(c("x","y","value","z"))

#' Coerce a NeuroSlice (or matrix-like) to a numeric matrix
#' @keywords internal
#' @noRd
slice_to_matrix <- function(slc) {
  plain_matrix <- function(x) {
    dx <- dim(x)
    if (length(dx) != 2L) {
      stop("Expected a 2D object when coercing to a plain matrix.", call. = FALSE)
    }
    matrix(
      as.numeric(x),
      nrow = dx[1],
      ncol = dx[2],
      dimnames = dimnames(x)
    )
  }

  if (is.matrix(slc) && !isS4(slc)) {
    return(plain_matrix(slc))
  }
  m <- try(as.matrix(slc), silent = TRUE)
  if (!inherits(m, "try-error") && is.matrix(m)) return(plain_matrix(m))
  a <- try(as.array(slc), silent = TRUE)
  if (!inherits(a, "try-error")) {
    if (length(dim(a)) == 2) return(plain_matrix(a))
    if (length(dim(a)) == 3) return(plain_matrix(a[, , 1, drop = TRUE]))
  }
  # try common slots / accessors without committing to class internals
  data <- try(slc@data, silent = TRUE)
  if (!inherits(data, "try-error") && is.matrix(data)) return(plain_matrix(data))
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

#' Orient slice-aligned matrices for raster rendering in world coordinates
#' @param slc A NeuroSlice describing the 2D slice geometry.
#' @param mat Numeric matrix aligned with \code{slc}.
#' @param alpha_map Optional numeric matrix aligned with \code{slc}.
#' @keywords internal
#' @noRd
orient_slice_for_raster <- function(slc, mat, alpha_map = NULL) {
  mat <- slice_to_matrix(mat)
  if (!is.null(alpha_map)) {
    alpha_map <- slice_to_matrix(alpha_map)
    if (!identical(dim(alpha_map), dim(mat))) {
      stop("'alpha_map' must have the same dimensions as 'mat'.", call. = FALSE)
    }
  }

  coords <- index_to_coord(space(slc), seq_len(length(slc)))
  xvals <- sort(unique(coords[, 1]))
  yvals <- sort(unique(coords[, 2]), decreasing = TRUE)

  out_mat <- matrix(NA_real_, nrow = length(yvals), ncol = length(xvals))
  out_alpha <- if (is.null(alpha_map)) NULL else {
    matrix(NA_real_, nrow = length(yvals), ncol = length(xvals))
  }

  col_idx <- match(coords[, 1], xvals)
  row_idx <- match(coords[, 2], yvals)
  fill_idx <- cbind(row_idx, col_idx)

  out_mat[fill_idx] <- as.numeric(mat)
  if (!is.null(out_alpha)) {
    out_alpha[fill_idx] <- as.numeric(alpha_map)
  }

  list(
    mat = out_mat,
    alpha_map = out_alpha,
    x = xvals,
    y = yvals
  )
}

#' Compute raster annotation extents from pixel centers
#' @param centers Numeric vector of pixel centers.
#' @keywords internal
#' @noRd
raster_extent_from_centers <- function(centers) {
  vals <- sort(unique(centers))
  if (!length(vals)) return(c(0, 1))
  if (length(vals) == 1L) return(c(vals[1] - 0.5, vals[1] + 0.5))
  step <- stats::median(diff(vals))
  c(vals[1] - step / 2, vals[length(vals)] + step / 2)
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
#' @param mat Numeric matrix.
#' @param cmap Palette name.
#' @param limits Numeric length-2 vector of display limits.
#' @param alpha Global alpha scalar (0..1).
#' @param alpha_map Optional numeric matrix (same dims as \code{mat}, values in
#'   \code{[0,1]}) providing per-pixel alpha. Effective alpha per pixel is
#'   \code{alpha_map[i] * alpha}. When \code{NULL} the scalar \code{alpha} is
#'   applied to all pixels.
#' @keywords internal
#' @noRd
matrix_to_colors <- function(mat, cmap = "grays", limits = NULL, alpha = 1,
                              alpha_map = NULL) {
  cols <- resolve_cmap(cmap, 256)
  if (is.null(limits)) limits <- range(mat[is.finite(mat)])
  s <- rescale01(as.numeric(mat), limits)
  s[!is.finite(s)] <- NA_real_
  # index into palette
  idx <- 1 + floor(s * (length(cols) - 1))
  col <- cols[pmax(1, pmin(length(cols), idx))]
  # apply alpha
  if (is.null(alpha_map)) {
    grDevices::adjustcolor(col, alpha.f = alpha)
  } else {
    eff <- as.numeric(alpha_map) * alpha
    mapply(function(c, a) grDevices::adjustcolor(c, alpha.f = a),
           col, eff, USE.NAMES = FALSE)
  }
}

#' Build a numeric RGBA array for a matrix
#' @param mat Numeric matrix.
#' @param cmap Palette name.
#' @param limits Numeric length-2 display limits.
#' @param alpha Global alpha scalar (0..1).
#' @param alpha_map Optional per-pixel alpha matrix (same dims as \code{mat},
#'   values in \code{[0,1]}). Effective alpha is \code{alpha_map * alpha}.
#' @keywords internal
#' @noRd
matrix_to_rgba <- function(mat, cmap = "grays", limits = NULL, alpha = 1,
                           alpha_map = NULL) {
  cols <- resolve_cmap(cmap, 256)
  if (is.null(limits)) limits <- range(mat[is.finite(mat)])

  s <- rescale01(as.numeric(mat), limits)
  idx <- 1 + floor(s * (length(cols) - 1))
  idx[!is.finite(idx)] <- 1L
  idx <- pmax(1L, pmin(length(cols), idx))

  rgb <- grDevices::col2rgb(cols[idx], alpha = TRUE) / 255

  if (is.null(alpha_map)) {
    eff_alpha <- rep(alpha, length(s))
  } else {
    eff_alpha <- as.numeric(alpha_map) * alpha
  }
  eff_alpha[!is.finite(eff_alpha)] <- 0
  eff_alpha <- pmax(0, pmin(1, eff_alpha))
  eff_alpha[!is.finite(s)] <- 0

  rgba <- array(0, dim = c(nrow(mat), ncol(mat), 4L))
  rgba[, , 1] <- matrix(rgb[1, ], nrow = nrow(mat), ncol = ncol(mat))
  rgba[, , 2] <- matrix(rgb[2, ], nrow = nrow(mat), ncol = ncol(mat))
  rgba[, , 3] <- matrix(rgb[3, ], nrow = nrow(mat), ncol = ncol(mat))
  rgba[, , 4] <- matrix(eff_alpha, nrow = nrow(mat), ncol = ncol(mat))

  rgba
}

#' Create a rasterGrob from a numeric matrix using a palette
#' @param mat Numeric matrix.
#' @param cmap Palette name.
#' @param limits Numeric length-2 display limits.
#' @param alpha Global alpha scalar (0..1).
#' @param alpha_map Optional per-pixel alpha matrix (same dims as \code{mat},
#'   values in \code{[0,1]}). Effective alpha is \code{alpha_map * alpha}.
#' @keywords internal
#' @noRd
matrix_to_raster_grob <- function(mat, cmap = "grays", limits = NULL, alpha = 1,
                                   alpha_map = NULL) {
  rgba <- matrix_to_rgba(
    mat = mat,
    cmap = cmap,
    limits = limits,
    alpha = alpha,
    alpha_map = alpha_map
  )
  grid::rasterGrob(rgba, interpolate = FALSE)
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
