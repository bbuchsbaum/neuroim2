#' Space utility functions
#'
#' Utilities for reasoning about mapped voxel spaces and slice embeddings.
#'
#' Compared to NiBabel-style helpers, these functions add a few R-friendly
#' improvements:
#' \itemize{
#'   \item Accept `NeuroSpace`, `NeuroVol`, and list/object `(shape, affine)` inputs.
#'   \item Handle inputs with more than 3 dimensions by using first 3 spatial dims.
#'   \item Support both R-style (1-based) and zero-based slice indexing.
#' }
#'
#' @name space_utils
#' @examples
#' sp <- NeuroSpace(c(10L, 8L, 6L), spacing = c(2, 2, 2))
#' out <- output_aligned_space(sp)
#' out$shape
#' out$affine
#'
#' slice_aff <- slice_to_volume_affine(index = 3, axis = 3, shape = c(10, 8, 6))
#' slice_aff
NULL

#' Parse mapped voxel input as `(shape, affine)`
#'
#' @param mapped_voxels Input object.
#' @return List with `shape` and `affine`.
#' @keywords internal
#' @noRd
.parse_mapped_voxels <- function(mapped_voxels) {
  if (inherits(mapped_voxels, "NeuroSpace")) {
    return(list(shape = dim(mapped_voxels), affine = trans(mapped_voxels)))
  }

  if (inherits(mapped_voxels, "NeuroVol") || inherits(mapped_voxels, "NeuroVec")) {
    return(list(shape = dim(mapped_voxels), affine = trans(space(mapped_voxels))))
  }

  if (is.list(mapped_voxels) && length(mapped_voxels) == 2) {
    if (!is.null(mapped_voxels$shape) && !is.null(mapped_voxels$affine)) {
      return(list(shape = mapped_voxels$shape, affine = mapped_voxels$affine))
    }
    return(list(shape = mapped_voxels[[1]], affine = mapped_voxels[[2]]))
  }

  if (!is.null(mapped_voxels$shape) && !is.null(mapped_voxels$affine)) {
    return(list(shape = mapped_voxels$shape, affine = mapped_voxels$affine))
  }

  stop("`mapped_voxels` must be a NeuroSpace/NeuroVol/NeuroVec or provide shape+affine")
}

#' Coerce affine to 4x4
#'
#' @param affine Affine matrix-like.
#' @param n_axes Number of spatial axes represented.
#' @return 4x4 affine matrix.
#' @keywords internal
#' @noRd
.as_affine4 <- function(affine, n_axes) {
  aff <- as.matrix(affine)
  if (!is.numeric(aff)) {
    stop("`affine` must be numeric")
  }

  if (all(dim(aff) == c(4L, 4L))) {
    return(aff)
  }

  expected <- n_axes + 1L
  if (!all(dim(aff) == c(expected, expected))) {
    stop("`affine` must be 4x4 or (n_axes+1)x(n_axes+1)")
  }

  out <- diag(4)
  out[seq_len(n_axes), seq_len(n_axes)] <- aff[seq_len(n_axes), seq_len(n_axes), drop = FALSE]
  out[seq_len(n_axes), 4] <- aff[seq_len(n_axes), expected]
  out
}

#' Build output-aligned voxel space for mapped voxels
#'
#' Returns a voxel space with positive diagonal affine that encloses all mapped
#' input voxel centers.
#'
#' @param mapped_voxels A `NeuroSpace`, `NeuroVol`, `NeuroVec`, object with
#'   `shape` and `affine`, or length-2 sequence `(shape, affine)`.
#' @param voxel_sizes Optional output voxel sizes for spatial axes.
#'   If scalar, treated as isotropic.
#' @return A named list with:
#'   \itemize{
#'     \item `shape`: output spatial shape.
#'     \item `affine`: output 4x4 affine (positive diagonal).
#'     \item `bounds`: world-space min/max of mapped corners.
#'   }
#' @rdname space_utils
#' @export
output_aligned_space <- function(mapped_voxels, voxel_sizes = NULL) {
  parsed <- .parse_mapped_voxels(mapped_voxels)

  shape <- as.numeric(parsed$shape)
  if (!is.numeric(shape) || length(shape) < 1) {
    stop("`shape` must be a numeric vector with at least one dimension")
  }
  if (any(!is.finite(shape)) || any(shape <= 0) || any(shape %% 1 != 0)) {
    stop("`shape` must contain positive integers")
  }

  n_axes_full <- length(shape)
  if (n_axes_full > 3) {
    warning("Using first three dimensions as spatial axes; trailing dimensions are ignored")
  }
  n_axes <- min(3L, n_axes_full)
  shape_spatial <- as.integer(shape[seq_len(n_axes)])

  aff4 <- .as_affine4(parsed$affine, n_axes = n_axes)

  out_vox <- rep(1, 3)
  if (!is.null(voxel_sizes)) {
    voxel_sizes <- as.numeric(voxel_sizes)
    if (length(voxel_sizes) == 1L) {
      voxel_sizes <- rep(voxel_sizes, n_axes)
    }
    if (length(voxel_sizes) != n_axes) {
      stop("`voxel_sizes` must have length 1 or match the number of spatial axes")
    }
    if (any(!is.finite(voxel_sizes)) || any(voxel_sizes <= 0)) {
      stop("`voxel_sizes` must contain strictly positive finite values")
    }
    out_vox[seq_len(n_axes)] <- voxel_sizes
  }

  # Corner coordinates in zero-based voxel index space.
  corner_ranges <- lapply(shape_spatial, function(n) c(0, n - 1))
  corners <- do.call(expand.grid, corner_ranges)
  corners <- as.matrix(corners)
  if (n_axes < 3L) {
    corners <- cbind(corners, matrix(0, nrow(corners), 3L - n_axes))
  }

  out_corners <- cbind(corners, 1) %*% t(aff4)
  out_corners <- out_corners[, 1:3, drop = FALSE]

  out_min <- apply(out_corners, 2, min)
  out_max <- apply(out_corners, 2, max)

  out_shape_full <- ceiling((out_max - out_min) / out_vox) + 1
  out_aff <- diag(c(out_vox, 1))
  out_aff[1:3, 4] <- out_min

  list(
    shape = as.integer(out_shape_full[seq_len(n_axes)]),
    affine = out_aff,
    bounds = list(min = out_min, max = out_max)
  )
}

#' Compatibility alias for `output_aligned_space()`
#'
#' @return Same as `output_aligned_space()`.
#' @rdname space_utils
#' @export
vox2out_vox <- function(mapped_voxels, voxel_sizes = NULL) {
  output_aligned_space(mapped_voxels, voxel_sizes = voxel_sizes)
}

#' Affine embedding of a 2D slice within a 3D volume
#'
#' Returns the affine mapping slice coordinates to 3D volume coordinates.
#'
#' @param index Slice index.
#' @param axis Slice axis (`1..3` or `0..2`).
#' @param shape Optional volume shape for bounds validation.
#' @param index_base Either `"R"` (1-based, default) or `"zero"`.
#' @return A `4 x 3` affine matrix from slice coordinates to volume coordinates.
#' @rdname space_utils
#' @export
slice_to_volume_affine <- function(index, axis, shape = NULL, index_base = c("R", "zero")) {
  index_base <- match.arg(index_base)

  if (!is.numeric(index) || length(index) != 1L || !is.finite(index) || index %% 1 != 0) {
    stop("`index` must be a finite integer scalar")
  }
  index <- as.integer(index)

  if (!is.numeric(axis) || length(axis) != 1L || !is.finite(axis) || axis %% 1 != 0) {
    stop("`axis` must be an integer scalar")
  }
  axis <- as.integer(axis)
  if (axis %in% 0:2) {
    axis <- axis + 1L
  }
  if (!(axis %in% 1:3)) {
    stop("`axis` must be in 1..3 (or 0..2)")
  }

  if (index_base == "R") {
    if (index < 1L) {
      stop("For `index_base = 'R'`, `index` must be >= 1")
    }
    index0 <- index - 1L
  } else {
    if (index < 0L) {
      stop("For `index_base = 'zero'`, `index` must be >= 0")
    }
    index0 <- index
  }

  if (!is.null(shape)) {
    shape <- as.numeric(shape)
    if (length(shape) < 3L || any(shape[1:3] <= 0) || any(shape[1:3] %% 1 != 0)) {
      stop("`shape` must provide at least 3 positive integer dimensions")
    }
    axis_n <- as.integer(shape[axis])
    if (index_base == "R" && index > axis_n) {
      stop("`index` exceeds the provided `shape` along `axis`")
    }
    if (index_base == "zero" && index >= axis_n) {
      stop("`index` exceeds the provided `shape` along `axis`")
    }
  }

  axes <- 1:4
  axes <- axes[axes != axis]
  slice_aff <- diag(4)[, axes, drop = FALSE]
  slice_aff[axis, 3] <- index0
  slice_aff
}

#' Compatibility alias for `slice_to_volume_affine()`
#'
#' @return Same as `slice_to_volume_affine()`.
#' @rdname space_utils
#' @export
slice2volume <- function(index, axis, shape = NULL, index_base = c("R", "zero")) {
  slice_to_volume_affine(index = index, axis = axis, shape = shape, index_base = index_base)
}
