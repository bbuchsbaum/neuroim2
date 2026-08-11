#' @include all_class.R
NULL
#' @include all_generic.R
NULL

#' Get Orientation Axis Codes
#'
#' Returns a character vector of anatomical axis direction labels for a
#' neuroimaging object or space. For example, \code{c("R", "A", "S")} for a
#' standard RAS-oriented image.
#'
#' @param x A \code{\linkS4class{NeuroVol}}, \code{\linkS4class{NeuroVec}},
#'   \code{\linkS4class{NeuroSpace}}, or a 4x4 affine matrix.
#' @return A character vector of length 3 with axis direction codes.
#'
#' @examples
#' sp <- NeuroSpace(c(10L, 10L, 10L), c(2, 2, 2))
#' axcodes(sp)
#'
#' vol <- DenseNeuroVol(array(0, c(10,10,10)), sp)
#' axcodes(vol)
#'
#' @rdname axcodes-methods
#' @export
setGeneric("axcodes", function(x) standardGeneric("axcodes"))

#' @rdname axcodes-methods
#' @export
setMethod("axcodes", "NeuroSpace", function(x) {
  affine_to_axcodes(trans(x))
})

#' @rdname axcodes-methods
#' @export
setMethod("axcodes", "NeuroObj", function(x) {
  affine_to_axcodes(trans(space(x)))
})

#' @rdname axcodes-methods
#' @export
setMethod("axcodes", "matrix", function(x) {
  affine_to_axcodes(x)
})


#' Reorient Image to Canonical (RAS+) Orientation
#'
#' Reorients a neuroimaging volume or vector to the canonical RAS+
#' (Right-Anterior-Superior) orientation by permuting and flipping axes.
#' This is equivalent to nibabel's \code{as_closest_canonical()}.
#'
#' @param x A \code{\linkS4class{NeuroVol}} or \code{\linkS4class{NeuroVec}}
#'   object.
#' @param target Character vector of length 3 giving the desired orientation.
#'   Default is \code{c("R", "A", "S")} (RAS+).
#' @return A reoriented object of the same class as \code{x}.
#'
#' @details
#' Reorienting to a target that differs only in which anatomical direction each
#' axis runs is a relabelling of the grid, not a resampling: the voxels are the
#' same voxels, permuted and flipped. This function does exactly that, so it is
#' exact, works on images of any size, and costs a single \code{aperm()}.
#'
#' It used to build the target space and hand it to \code{\link{resample}},
#' which had three consequences: values were interpolated where they should
#' merely have moved, the registration backend refused images smaller than four
#' voxels on any axis, and -- because the target space kept the source's
#' dimensions -- an axis-permuting reorientation wrote into a grid that did not
#' contain the data and silently lost about a quarter of it.
#'
#' @examples
#' sp <- NeuroSpace(c(10L, 10L, 10L), c(2, 2, 2))
#' vol <- DenseNeuroVol(array(rnorm(1000), c(10,10,10)), sp)
#' ras_vol <- as_canonical(vol)
#' axcodes(ras_vol)  # "R" "A" "S"
#'
#' @seealso \code{\link{axcodes}}, \code{\link{reorient}},
#'   \code{\link{apply_orientation}}
#' @export
as_canonical <- function(x, target = c("R", "A", "S")) {
  stopifnot(length(target) == 3)
  if (!inherits(x, "NeuroVol") && !inherits(x, "NeuroVec")) {
    cli::cli_abort("{.fn as_canonical} requires a {.cls NeuroVol} or {.cls NeuroVec} object.")
  }
  current <- axcodes(x)
  if (identical(current, target)) {
    return(x)
  }

  ornt <- .reorient_transform(space(x), target)
  # apply_orientation() permutes and flips the leading three axes and leaves
  # any trailing one (time) alone, so the same call serves both classes.
  arr <- apply_orientation(as.array(x), ornt)

  if (inherits(x, "NeuroVol")) {
    DenseNeuroVol(arr, reorient(space(x), target))
  } else {
    d <- dim(x)
    sp3 <- reorient(drop_dim(space(x)), target)
    DenseNeuroVec(arr, add_dim(sp3, d[4]),
                  label = x@label, volume_labels = volume_labels(x))
  }
}
