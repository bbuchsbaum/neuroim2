#' Affine utility functions
#'
#' Utilities for point transforms and affine decomposition/composition.
#'
#' These functions mirror core NiBabel affine helpers while using R-first
#' conventions and stronger argument checks.
#'
#' @name affine_utils
#' @examples
#' aff <- diag(c(2, 3, 4, 1))
#' aff[1:3, 4] <- c(10, 20, 30)
#'
#' pts <- rbind(c(1, 2, 3), c(4, 5, 6))
#' apply_affine(aff, pts)
#'
#' mv <- to_matvec(aff)
#' from_matvec(mv$matrix, mv$vector)
NULL

#' Raise affine error with a specific condition class
#'
#' @param message Error text.
#' @keywords internal
#' @noRd
.affine_error <- function(message) {
  err <- structure(list(message = message, call = NULL), class = c("AffineError", "error", "condition"))
  stop(err)
}

#' Validate homogeneous affine matrix
#'
#' @param aff Matrix-like affine transform.
#' @return Numeric matrix.
#' @keywords internal
#' @noRd
.validate_affine <- function(aff) {
  aff <- as.matrix(aff)
  if (!is.numeric(aff) || length(dim(aff)) != 2 || nrow(aff) < 2 || ncol(aff) < 2) {
    .affine_error("`aff` must be a numeric matrix with at least 2 rows and 2 columns")
  }
  aff
}

#' Apply an affine transform to points
#'
#' Applies a homogeneous affine transform to points where the coordinate axis is
#' the last dimension.
#'
#' @param aff `(N x N)` affine matrix in homogeneous coordinates.
#' @param pts Points with last dimension `N - 1`; shape may be vector, matrix,
#'   or higher-dimensional array.
#' @param inplace Included for API compatibility; ignored in R.
#' @return Points transformed by `aff`, with same leading shape as `pts`.
#' @rdname affine_utils
#' @export
apply_affine <- function(aff, pts, inplace = FALSE) {
  aff <- .validate_affine(aff)
  pts <- as.array(pts)
  dims <- dim(pts)

  nd_in <- ncol(aff) - 1L
  nd_out <- nrow(aff) - 1L

  if (is.null(dims)) {
    dims <- length(pts)
    pts <- array(pts, dim = dims)
  }

  if (length(dims) == 1L) {
    if (dims[1] != nd_in) {
      .affine_error("For vector `pts`, length must match affine input dimension")
    }
    pts_mat <- matrix(pts, nrow = 1L, ncol = nd_in)
    leading_dims <- integer(0)
  } else {
    if (dims[length(dims)] != nd_in) {
      .affine_error("Last dimension of `pts` must match affine input dimension")
    }
    leading_dims <- dims[-length(dims)]
    perm <- c(length(dims), seq_len(length(dims) - 1L))
    pts_perm <- aperm(pts, perm)
    n_pts <- prod(leading_dims)
    pts_mat <- t(matrix(pts_perm, nrow = nd_in, ncol = n_pts))
  }

  rzs <- aff[seq_len(nd_out), seq_len(nd_in), drop = FALSE]
  trans <- aff[seq_len(nd_out), ncol(aff), drop = TRUE]
  out_mat <- pts_mat %*% t(rzs) + matrix(trans, nrow = nrow(pts_mat), ncol = nd_out, byrow = TRUE)

  if (length(dims) == 1L) {
    return(as.numeric(out_mat[1, ]))
  }

  out_arr <- array(t(out_mat), dim = c(nd_out, leading_dims))
  if (length(leading_dims) == 1L) {
    return(t(out_arr))
  }
  aperm(out_arr, c(seq_len(length(leading_dims)) + 1L, 1L))
}

#' Split affine into linear matrix and translation vector
#'
#' @param transform Homogeneous transform matrix.
#' @return A list with `matrix` and `vector`.
#' @rdname affine_utils
#' @export
to_matvec <- function(transform) {
  transform <- .validate_affine(transform)
  nd_out <- nrow(transform) - 1L
  nd_in <- ncol(transform) - 1L

  list(
    matrix = transform[seq_len(nd_out), seq_len(nd_in), drop = FALSE],
    vector = transform[seq_len(nd_out), ncol(transform), drop = TRUE]
  )
}

#' Combine matrix and translation vector into a homogeneous affine
#'
#' @param matrix Linear part of transform (`N x M`).
#' @param vector Optional translation vector of length `N`.
#' @return Homogeneous affine matrix of shape `(N+1) x (M+1)`.
#' @rdname affine_utils
#' @export
from_matvec <- function(matrix, vector = NULL) {
  matrix <- as.matrix(matrix)
  if (!is.numeric(matrix) || length(dim(matrix)) != 2) {
    .affine_error("`matrix` must be a numeric 2D matrix")
  }

  n_out <- nrow(matrix)
  n_in <- ncol(matrix)
  out <- matrix(0, nrow = n_out + 1L, ncol = n_in + 1L)
  out[seq_len(n_out), seq_len(n_in)] <- matrix
  out[n_out + 1L, n_in + 1L] <- 1

  if (!is.null(vector)) {
    vector <- as.numeric(vector)
    if (length(vector) != n_out || any(!is.finite(vector))) {
      .affine_error("`vector` must be finite numeric with length equal to nrow(matrix)")
    }
    out[seq_len(n_out), n_in + 1L] <- vector
  }

  out
}

#' Append diagonal axes to an affine
#'
#' Useful when expanding an affine to include additional dimensions.
#'
#' @param aff Base affine matrix.
#' @param steps Diagonal values for appended dimensions.
#' @param starts Optional translations for appended dimensions.
#' @return Expanded affine matrix.
#' @rdname affine_utils
#' @export
append_diag <- function(aff, steps, starts = numeric(0)) {
  aff <- .validate_affine(aff)
  steps <- as.numeric(steps)
  starts <- as.numeric(starts)

  if (length(steps) < 1L || any(!is.finite(steps))) {
    .affine_error("`steps` must contain at least one finite numeric value")
  }

  n_steps <- length(steps)
  if (length(starts) == 0L) {
    starts <- rep(0, n_steps)
  } else if (length(starts) != n_steps || any(!is.finite(starts))) {
    .affine_error("`starts` must be empty or have same finite length as `steps`")
  }

  old_n_out <- nrow(aff) - 1L
  old_n_in <- ncol(aff) - 1L

  out <- matrix(
    0,
    nrow = old_n_out + n_steps + 1L,
    ncol = old_n_in + n_steps + 1L
  )

  out[seq_len(old_n_out), seq_len(old_n_in)] <- aff[seq_len(old_n_out), seq_len(old_n_in), drop = FALSE]
  out[seq_len(old_n_out), ncol(out)] <- aff[seq_len(old_n_out), ncol(aff), drop = TRUE]

  for (i in seq_len(n_steps)) {
    out[old_n_out + i, old_n_in + i] <- steps[i]
  }

  out[(old_n_out + 1L):(old_n_out + n_steps), ncol(out)] <- starts
  out[nrow(out), ncol(out)] <- 1
  out
}

#' Right-associated matrix product reduction
#'
#' For arrays `A, B, C, ...`, returns `A %*% (B %*% (C %*% ...))`.
#'
#' @param ... Matrices/vectors compatible with \code{\%*\%}.
#' @return Matrix product.
#' @rdname affine_utils
#' @export
dot_reduce <- function(...) {
  mats <- list(...)
  if (length(mats) < 1L) {
    .affine_error("`dot_reduce()` needs at least one argument")
  }
  Reduce(function(x, y) y %*% x, rev(mats))
}

#' Compute voxel sizes from affine columns
#'
#' @param affine Affine matrix.
#' @return Numeric vector of voxel sizes (column norms of linear block).
#' @rdname affine_utils
#' @export
voxel_sizes <- function(affine) {
  affine <- .validate_affine(affine)
  top_left <- affine[-nrow(affine), -ncol(affine), drop = FALSE]
  sqrt(colSums(top_left * top_left))
}

#' Estimate affine obliquity angles
#'
#' Returns per-output-axis obliquity relative to cardinal axes, in radians.
#'
#' @param affine Affine matrix.
#' @return Numeric vector of obliquity angles.
#' @rdname affine_utils
#' @export
obliquity <- function(affine) {
  affine <- .validate_affine(affine)
  rzs <- affine[-nrow(affine), -ncol(affine), drop = FALSE]
  vs <- voxel_sizes(affine)

  if (any(vs == 0)) {
    .affine_error("Cannot compute obliquity for zero-length voxel axes")
  }

  cos_mat <- abs(sweep(rzs, 2, vs, "/"))
  best_cosines <- apply(cos_mat, 1, max)
  acos(pmax(pmin(best_cosines, 1), -1))
}

#' Rescale affine to new voxel sizes while preserving image center
#'
#' Preserves rotations/shears and world location of the center voxel.
#'
#' @param affine Square homogeneous affine matrix.
#' @param shape Original grid shape (spatial axes).
#' @param zooms New voxel sizes for spatial axes.
#' @param new_shape Optional new grid shape; defaults to `shape`.
#' @return Rescaled affine matrix.
#' @rdname affine_utils
#' @export
rescale_affine <- function(affine, shape, zooms, new_shape = NULL) {
  affine <- .validate_affine(affine)
  if (nrow(affine) != ncol(affine)) {
    .affine_error("`rescale_affine()` requires a square homogeneous affine")
  }

  ndim <- nrow(affine) - 1L
  shape <- as.numeric(shape)
  zooms <- as.numeric(zooms)
  new_shape <- as.numeric(if (is.null(new_shape)) shape else new_shape)

  if (length(shape) != ndim || length(zooms) != ndim || length(new_shape) != ndim) {
    .affine_error("`shape`, `zooms`, and `new_shape` must all match affine spatial dimensionality")
  }
  if (any(shape <= 0) || any(new_shape <= 0) || any(zooms <= 0)) {
    .affine_error("`shape`, `new_shape`, and `zooms` must be strictly positive")
  }

  s <- voxel_sizes(affine)
  if (any(s == 0)) {
    .affine_error("Cannot rescale affine with zero voxel size in linear block")
  }

  rzs <- affine[seq_len(ndim), seq_len(ndim), drop = FALSE]
  rzs_out <- sweep(rzs, 2, zooms / s, "*")

  # Match NiBabel semantics: use integer center voxel via floor.
  center_in <- floor((shape - 1) / 2)
  center_out <- floor((new_shape - 1) / 2)

  centroid <- as.numeric(apply_affine(affine, center_in))
  t_out <- centroid - as.numeric(rzs_out %*% center_out)

  from_matvec(rzs_out, t_out)
}
