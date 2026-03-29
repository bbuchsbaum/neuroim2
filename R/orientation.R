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
#' The function works by computing the orientation transform from the
#' current axis codes to the target codes, then applying the necessary
#' axis permutations and flips. The affine matrix is updated to reflect
#' the new orientation while preserving world-coordinate mapping.
#'
#' @examples
#' sp <- NeuroSpace(c(10L, 10L, 10L), c(2, 2, 2))
#' vol <- DenseNeuroVol(array(rnorm(1000), c(10,10,10)), sp)
#' ras_vol <- as_canonical(vol)
#' axcodes(ras_vol)  # "R" "A" "S"
#'
#' @seealso \code{\link{axcodes}}, \code{\link{reorient}}
#' @export
as_canonical <- function(x, target = c("R", "A", "S")) {
  stopifnot(length(target) == 3)
  current <- axcodes(x)
  if (identical(current, target)) {
    return(x)
  }

  # Use existing reorient infrastructure
  if (inherits(x, "NeuroVol")) {
    new_space <- reorient(space(x), target)
    resample(x, new_space)
  } else if (inherits(x, "NeuroVec")) {
    new_space <- reorient(space(x), target)
    # Resample each volume and reconstruct
    d <- dim(x)
    vols <- lapply(seq_len(d[4]), function(i) {
      resample(x[[i]], new_space)
    })
    mat <- do.call(cbind, lapply(vols, function(v) as.numeric(v@.Data)))
    dspace <- add_dim(new_space, d[4])
    DenseNeuroVec(mat, dspace)
  } else {
    cli::cli_abort("{.fn as_canonical} requires a {.cls NeuroVol} or {.cls NeuroVec} object.")
  }
}
