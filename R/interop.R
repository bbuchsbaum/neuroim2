#' @include all_class.R all_generic.R
NULL

#' Convert between neuroim2 objects and RNifti / oro.nifti images
#'
#' @description
#' \code{as_neuro()} turns an image held by another R neuroimaging package
#' into a neuroim2 object, and \code{as_nifti()} goes the other way. Together
#' they let neuroim2 take part in a pipeline built on \pkg{RNifti},
#' \pkg{RNiftyReg}, \pkg{oro.nifti}, \pkg{neurobase} or \pkg{fslr} without a
#' round trip through a file.
#'
#' @details
#' \code{as_neuro()} accepts
#' \itemize{
#'   \item an \pkg{RNifti} \code{niftiImage} or \code{internalImage}, e.g. the
#'     result of \code{RNifti::readNifti()} or of an \pkg{RNiftyReg}
#'     registration;
#'   \item an \pkg{oro.nifti} \code{nifti} object, e.g. the result of
#'     \code{oro.nifti::readNIfTI()}, \code{neurobase::readnii()} or an
#'     \pkg{fslr} command;
#'   \item a file name, which is passed to \code{\link{read_image}};
#'   \item a neuroim2 object, which is returned unchanged.
#' }
#' A 3-D image becomes a \code{\linkS4class{DenseNeuroVol}} and a 4-D image
#' a \code{\linkS4class{DenseNeuroVec}}. Trailing dimensions of extent 1 are
#' dropped first, so a 4-D image with one volume becomes a volume.
#'
#' The voxel-to-world affine is taken from \code{RNifti::xform()}, which
#' applies the NIfTI rule: the sform if it is set, else the qform, else a
#' scaling by the voxel sizes. \pkg{RNifti} and neuroim2 share one convention
#' for it (zero-based voxel indices to RAS+ millimetres), so
#' \code{trans(as_neuro(img))} equals \code{RNifti::xform(img)} to the seven
#' significant figures a \code{\linkS4class{NeuroSpace}} keeps. Data scaling
#' (\code{scl_slope}, \code{scl_inter}) has already been applied by the
#' package that read the image, so values carry over as they are.
#'
#' \code{as_nifti()} returns an \pkg{RNifti} \code{niftiImage} holding the
#' values of \code{x} as a dense array. The affine of \code{x} is written to
#' both the sform and the qform, and the voxel sizes to \code{pixdim}. The
#' qform can only describe a rotation with scaling, so for an affine with
#' shear the qform is the closest such transform and the sform is the exact
#' one; readers that follow the NIfTI rule use the sform. Any object with an
#' \code{as.array()} method and a \code{\link{space}} can be converted, so
#' sparse, clustered, memory-mapped and file-backed objects are read into
#' memory by the call. A logical volume is stored as 0/1 integers.
#'
#' Because \code{RNifti::asNifti()} is an S3 generic, \code{asNifti(x)} works
#' on neuroim2 objects as well and gives the same result as
#' \code{as_nifti(x)}. To get an \pkg{oro.nifti} object, pass the result on:
#' \code{oro.nifti::nii2oro(as_nifti(x))}.
#'
#' Volume labels, NIfTI extensions and header fields other than the geometry
#' are not carried across in either direction.
#'
#' @param x For \code{as_neuro()}, a \code{niftiImage}, \code{internalImage},
#'   \pkg{oro.nifti} \code{nifti} object, file name or neuroim2 object. For
#'   \code{as_nifti()}, a \code{\linkS4class{NeuroVol}} or
#'   \code{\linkS4class{NeuroVec}}.
#' @param ... For \code{as_neuro()}, further arguments passed to
#'   \code{\link{read_image}} when \code{x} is a file name; otherwise unused.
#' @param label For \code{as_neuro()}, an optional label for the result.
#' @param xform_code The code written to \code{sform_code} and
#'   \code{qform_code}: 1 (scanner anatomical, the default), 2 (aligned to
#'   another image), 3 (Talairach) or 4 (MNI152). neuroim2 does not record
#'   which of these a space is, so say so here when it matters downstream.
#' @param reference Ignored; present because \code{RNifti::asNifti()} has it.
#'
#' @return \code{as_neuro()} returns a \code{\linkS4class{DenseNeuroVol}} or
#'   \code{\linkS4class{DenseNeuroVec}}. \code{as_nifti()} returns a
#'   \code{niftiImage}.
#'
#' @seealso \code{\link{read_image}}, \code{\link{write_vol}},
#'   \code{\link{as_nifti_header}}
#'
#' @examples
#' fname <- system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2")
#'
#' # An image read by RNifti becomes a neuroim2 volume ...
#' img <- RNifti::readNifti(fname)
#' vol <- as_neuro(img)
#' all.equal(trans(vol), RNifti::xform(img), check.attributes = FALSE,
#'           tolerance = 1e-6)
#'
#' # ... and a neuroim2 volume goes back, with its geometry
#' smoothed <- gaussian_blur(vol, as.mask(vol > 0), sigma = 2)
#' img2 <- as_nifti(smoothed)
#' RNifti::pixdim(img2)
#'
#' @name as_neuro
#' @rdname as_neuro
#' @export
as_neuro <- function(x, ..., label = "") {
  if (is(x, "NeuroObj")) {
    return(x)
  }

  # An internalImage is a character vector underneath, so test for a file
  # name only once images are ruled out.
  if (is.character(x) && !inherits(x, "niftiImage")) {
    return(read_image(x, ...))
  }

  if (isS4(x) && is(x, "nifti")) {
    # An oro.nifti object; RNifti has the bridge.
    x <- RNifti::asNifti(x, internal = FALSE)
  }

  if (!inherits(x, "niftiImage")) {
    cli::cli_abort(c(
      "{.fn as_neuro} cannot convert an object of class {.cls {class(x)[1]}}.",
      "i" = "It takes an {.pkg RNifti} {.cls niftiImage}, an {.pkg oro.nifti} {.cls nifti} object or a file name.",
      "i" = "For a plain array, give it a space: {.code NeuroVol(x, NeuroSpace(dim(x), spacing))}."
    ))
  }

  aff <- unclass(RNifti::xform(x))
  attributes(aff) <- list(dim = c(4L, 4L))

  dat <- if (inherits(x, "internalImage")) as.array(x) else x
  d <- dim(dat)
  if (inherits(dat, "rgbArray") || is.complex(dat)) {
    cli::cli_abort("{.fn as_neuro} does not convert RGB or complex images.")
  }
  attributes(dat) <- list(dim = d)

  # Trailing singleton dimensions carry no information: NIfTI writers often
  # store a volume as x * y * z * 1.
  while (length(d) > 3L && d[length(d)] == 1L) {
    d <- d[-length(d)]
  }
  dim(dat) <- d

  nd <- length(d)
  if (nd < 3L || nd > 4L) {
    cli::cli_abort("{.fn as_neuro} converts 3-D and 4-D images, not {nd}-D.")
  }

  sp <- NeuroSpace(d, trans = aff)
  if (nd == 3L) {
    DenseNeuroVol(dat, sp, label = label)
  } else {
    DenseNeuroVec(dat, sp, label = if (nzchar(label)) label else "none")
  }
}

#' @rdname as_neuro
#' @export
as_nifti <- function(x, xform_code = 1L) {
  if (inherits(x, "niftiImage")) {
    return(x)
  }
  if (!is(x, "NeuroObj")) {
    cli::cli_abort(c(
      "{.fn as_nifti} converts neuroim2 volumes and vectors, not {.cls {class(x)[1]}}.",
      "i" = "Use {.fn RNifti::asNifti} for arrays and other image classes."
    ))
  }
  xform_code <- as.integer(xform_code)
  if (length(xform_code) != 1L || is.na(xform_code) || xform_code < 1L || xform_code > 4L) {
    cli::cli_abort("{.arg xform_code} must be 1, 2, 3 or 4.")
  }

  sp <- space(x)
  d <- dim(sp)
  if (length(d) < 3L || length(d) > 4L) {
    cli::cli_abort("{.fn as_nifti} converts 3-D and 4-D images, not {length(d)}-D.")
  }

  dat <- as.array(x)
  if (is.logical(dat)) {
    storage.mode(dat) <- "integer"
  }
  attributes(dat) <- list(dim = d)

  aff <- trans(sp)
  if (!identical(dim(aff), c(4L, 4L))) {
    cli::cli_abort("{.fn as_nifti} needs a 4 x 4 affine; this space has a {nrow(aff)} x {ncol(aff)} one.")
  }

  img <- RNifti::asNifti(dat)
  # pixdim first: the qform setter reads the voxel sizes back out of the
  # matrix, and the two must agree.
  pd <- rep(1, length(d))
  pd[1:3] <- sqrt(colSums(aff[1:3, 1:3]^2))
  RNifti::pixdim(img) <- pd
  RNifti::qform(img) <- structure(aff, code = xform_code)
  RNifti::sform(img) <- structure(aff, code = xform_code)
  img
}

#' @rdname as_neuro
#' @exportS3Method RNifti::asNifti
asNifti.NeuroObj <- function(x, reference = NULL, ..., xform_code = 1L) {
  as_nifti(x, xform_code = xform_code)
}
