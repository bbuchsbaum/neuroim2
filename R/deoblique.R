#' Deoblique a Neuroimaging Space or Volume
#'
#' AFNI-like helper that mirrors the core behavior of
#' \code{3dWarp -deoblique}:
#' \itemize{
#'   \item If \code{gridset} is supplied, use that as the output grid.
#'   \item Otherwise, build an axis-aligned output grid that encloses the
#'   input field-of-view.
#'   \item If \code{newgrid} is not supplied, use the minimum input voxel size
#'   isotropically (AFNI-style default for deobliquing).
#' }
#'
#' For \code{NeuroSpace}, this returns the target deobliqued space.
#' For \code{NeuroVol}, it also resamples image data into that space.
#'
#' @param x A \code{NeuroSpace} or \code{NeuroVol}.
#' @param gridset Optional output grid (a \code{NeuroSpace} or \code{NeuroVol}),
#'   analogous to AFNI's \code{-gridset}. Mutually exclusive with \code{newgrid}.
#' @param newgrid Optional scalar output voxel size (mm), analogous to AFNI's
#'   \code{-newgrid}. If omitted and \code{gridset} is \code{NULL}, the minimum
#'   input voxel size is used isotropically.
#' @param method Interpolation method used when \code{x} is a \code{NeuroVol}:
#'   one of \code{"nearest"}, \code{"linear"}, or \code{"cubic"}.
#' @param engine Resampling engine passed to \code{\link{resample_to}}.
#'
#' @return If \code{x} is a \code{NeuroSpace}, returns a deobliqued
#'   \code{NeuroSpace}. If \code{x} is a \code{NeuroVol}, returns a resampled
#'   \code{NeuroVol} in deobliqued space.
#'
#' @seealso \code{\link{output_aligned_space}}, \code{\link{resample_to}}
#'
#' @examples
#' sp <- NeuroSpace(c(32, 32, 20), spacing = c(2, 2, 3))
#' tx <- trans(sp)
#' tx[1, 2] <- 0.15
#' sp_obl <- NeuroSpace(dim(sp), spacing = spacing(sp), trans = tx)
#'
#' # Build deobliqued target space (minimum spacing default)
#' sp_deob <- deoblique(sp_obl)
#'
#' # Resample a volume to deobliqued space
#' vol <- NeuroVol(array(rnorm(prod(dim(sp_obl))), dim(sp_obl)), sp_obl)
#' \donttest{
#' vol_deob <- deoblique(vol, method = "linear")
#' }
#'
#' @export
deoblique <- function(x,
                      gridset = NULL,
                      newgrid = NULL,
                      method = c("linear", "nearest", "cubic"),
                      engine = c("internal")) {

  if (!inherits(x, "NeuroSpace") && !inherits(x, "NeuroVol")) {
    stop("`x` must be a NeuroSpace or NeuroVol.")
  }

  if (!is.null(gridset) && !is.null(newgrid)) {
    stop("`gridset` and `newgrid` are mutually exclusive.")
  }

  source_space <- if (inherits(x, "NeuroVol")) space(x) else x
  if (ndim(source_space) != 3L) {
    stop("`deoblique()` currently supports 3D spaces/volumes only.")
  }

  if (!is.null(gridset)) {
    target_space <- if (inherits(gridset, "NeuroVol")) space(gridset) else gridset
    if (!inherits(target_space, "NeuroSpace")) {
      stop("`gridset` must be a NeuroSpace or NeuroVol.")
    }
    if (ndim(target_space) != 3L) {
      stop("`gridset` must define a 3D NeuroSpace/NeuroVol.")
    }
  } else {
    if (is.null(newgrid)) {
      # AFNI-style deoblique default: isotropic grid at the minimum input spacing.
      newgrid <- min(spacing(source_space)[seq_len(min(3L, length(spacing(source_space))))])
    }

    if (!is.numeric(newgrid) || length(newgrid) != 1L || !is.finite(newgrid) || newgrid <= 0) {
      stop("`newgrid` must be a single positive numeric value.")
    }

    out <- output_aligned_space(source_space, voxel_sizes = as.numeric(newgrid))
    target_space <- NeuroSpace(
      dim = out$shape,
      spacing = diag(out$affine)[1:3],
      origin = out$affine[1:3, 4],
      trans = out$affine
    )
  }

  if (inherits(x, "NeuroSpace")) {
    return(target_space)
  }

  method <- match.arg(method)
  engine <- match.arg(engine, c("internal", "RNiftyReg"))
  resample_to(x, target_space, method = method, engine = engine)
}
