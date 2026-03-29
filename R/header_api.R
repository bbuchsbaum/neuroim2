#' @include all_class.R
NULL
#' @include all_generic.R
NULL

#' Access NIfTI Header Information
#'
#' @description
#' Retrieves header metadata from neuroimaging objects or files. Returns
#' a structured list with commonly needed fields like sform/qform matrices,
#' TR, intent codes, and data scaling parameters.
#'
#' For low-level access to all raw NIfTI header fields, use the
#' \code{$raw} element of the returned list.
#'
#' @param x A \code{\linkS4class{NeuroVol}}, \code{\linkS4class{NeuroVec}},
#'   \code{\linkS4class{FileMetaInfo}}, or a character file path.
#' @return A list of class \code{"NeuroHeader"} with elements:
#'   \describe{
#'     \item{dim}{Integer dimensions.}
#'     \item{pixdim}{Voxel sizes including TR.}
#'     \item{spacing}{Spatial voxel sizes (first 3 pixdim values).}
#'     \item{origin}{Coordinate origin.}
#'     \item{trans}{4x4 affine transform.}
#'     \item{qform}{List with \code{matrix} (4x4) and \code{code} (integer).}
#'     \item{sform}{List with \code{matrix} (4x4) and \code{code} (integer).}
#'     \item{intent_code}{NIfTI intent code.}
#'     \item{intent_name}{NIfTI intent name string.}
#'     \item{descrip}{Description string.}
#'     \item{data_type}{Storage data type label.}
#'     \item{bitpix}{Bits per pixel.}
#'     \item{scl_slope}{Data scaling slope.}
#'     \item{scl_inter}{Data scaling intercept.}
#'     \item{cal_min}{Display intensity minimum.}
#'     \item{cal_max}{Display intensity maximum.}
#'     \item{TR}{Repetition time (4th pixdim), or \code{NA} if not 4D.}
#'     \item{raw}{The complete raw header list (all NIfTI fields).}
#'   }
#'
#' @examples
#' \donttest{
#' f <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
#' h <- header(f)
#' h$dim
#' h$TR
#' h$sform
#' h$descrip
#' }
#'
#' @name header
#' @rdname header-methods
#' @export
setGeneric("header", function(x) standardGeneric("header"))

#' @keywords internal
#' @noRd
.build_neuro_header <- function(mi) {
  hdr <- if (.hasSlot(mi, "header") && is.list(mi@header)) mi@header else list()
  dims <- as.integer(dim(mi))

  # pixdim: from header if available, else reconstruct from spacing
  pixdim <- if (!is.null(hdr$pixdim)) as.numeric(hdr$pixdim) else {
    sp <- if (.hasSlot(mi, "spacing")) as.numeric(mi@spacing) else numeric(0)
    c(0, sp, rep(0, max(0, 4 - length(sp))))
  }

  spacing <- if (.hasSlot(mi, "spacing")) as.numeric(mi@spacing) else numeric(0)
  origin_val <- if (.hasSlot(mi, "origin")) as.numeric(mi@origin) else numeric(0)
  tx <- trans(mi)

  # sform: already a 4x4 matrix in the header list for NIFTIMetaInfo
  sform_code <- if (!is.null(hdr$sform_code)) as.integer(hdr$sform_code) else NA_integer_
  sform_mat <- if (is.matrix(hdr$sform) && nrow(hdr$sform) == 4 && ncol(hdr$sform) == 4) {
    hdr$sform
  } else {
    tx
  }

  # qform: already a 4x4 matrix in the header list for NIFTIMetaInfo
  qform_code <- if (!is.null(hdr$qform_code)) as.integer(hdr$qform_code) else NA_integer_
  qform_mat <- if (is.matrix(hdr$qform) && nrow(hdr$qform) == 4 && ncol(hdr$qform) == 4) {
    hdr$qform
  } else {
    NA
  }

  # TR: pixdim[5] is the 5th element (index 5), which is the time step for 4D
  tr_val <- if (length(dims) >= 4 && length(pixdim) >= 5) as.numeric(pixdim[5]) else NA_real_

  # Scaling: NIfTI uses scl_slope / scl_intercept in the header list
  slope <- if (!is.null(hdr$scl_slope)) as.numeric(hdr$scl_slope) else {
    if (.hasSlot(mi, "slope") && length(mi@slope) > 0) mi@slope[1] else 1
  }
  inter <- if (!is.null(hdr$scl_intercept)) as.numeric(hdr$scl_intercept) else {
    if (.hasSlot(mi, "intercept") && length(mi@intercept) > 0) mi@intercept[1] else 0
  }

  # Description: stored as $description in our header list.
  # The raw value may be a numeric byte vector (read as integers) or a
  # character vector; collapse and strip null terminators either way.
  descrip_val <- if (!is.null(hdr$description)) {
    d <- hdr$description
    if (is.numeric(d)) {
      # Treat as raw bytes: drop zeros and convert to string
      d <- d[d != 0]
      if (length(d) == 0) NA_character_ else rawToChar(as.raw(d))
    } else {
      s <- paste(as.character(d), collapse = "")
      if (nzchar(s)) s else NA_character_
    }
  } else if (!is.null(hdr$descrip)) {
    as.character(hdr$descrip)
  } else {
    NA_character_
  }

  # data_type label
  dtype_val <- if (.hasSlot(mi, "data_type") && length(mi@data_type) > 0) {
    mi@data_type
  } else if (!is.null(hdr$data_storage)) {
    as.character(hdr$data_storage)
  } else {
    NA_character_
  }

  # intent_name: stored as $intent_name (character vector of length 16 in raw)
  intent_name_val <- if (!is.null(hdr$intent_name)) {
    paste(as.character(hdr$intent_name), collapse = "")
  } else {
    NA_character_
  }

  out <- list(
    dim         = dims,
    pixdim      = pixdim,
    spacing     = spacing,
    origin      = origin_val,
    trans       = tx,
    qform       = list(matrix = qform_mat, code = qform_code),
    sform       = list(matrix = sform_mat, code = sform_code),
    intent_code = if (!is.null(hdr$intent_code)) as.integer(hdr$intent_code) else NA_integer_,
    intent_name = intent_name_val,
    descrip     = descrip_val,
    data_type   = dtype_val,
    bitpix      = if (!is.null(hdr$bitpix)) as.integer(hdr$bitpix) else NA_integer_,
    scl_slope   = slope,
    scl_inter   = inter,
    cal_min     = if (!is.null(hdr$cal_min)) as.numeric(hdr$cal_min) else NA_real_,
    cal_max     = if (!is.null(hdr$cal_max)) as.numeric(hdr$cal_max) else NA_real_,
    TR          = tr_val,
    raw         = hdr
  )
  class(out) <- "NeuroHeader"
  out
}

#' @rdname header-methods
#' @export
setMethod("header", signature(x = "FileMetaInfo"), function(x) {
  .build_neuro_header(x)
})

#' @rdname header-methods
#' @export
setMethod("header", signature(x = "character"), function(x) {
  mi <- read_header(x)
  .build_neuro_header(mi)
})

#' @method print NeuroHeader
#' @export
print.NeuroHeader <- function(x, ...) {
  cat(sprintf("NIfTI Header [%s]\n", paste(x$dim, collapse = " x ")))
  cat(sprintf("  Data type  : %s (%d-bit)\n", x$data_type,
              if (!is.na(x$bitpix)) x$bitpix else 0L))
  cat(sprintf("  Spacing    : %s mm\n",
              paste(round(x$spacing, 3), collapse = " x ")))
  cat(sprintf("  Origin     : %s\n",
              paste(round(x$origin, 2), collapse = ", ")))
  if (!is.na(x$sform$code)) cat(sprintf("  sform code : %d\n", x$sform$code))
  if (!is.na(x$qform$code)) cat(sprintf("  qform code : %d\n", x$qform$code))
  if (!is.na(x$TR) && x$TR > 0) cat(sprintf("  TR         : %.3f s\n", x$TR))
  if (!is.na(x$intent_code) && x$intent_code != 0L) {
    cat(sprintf("  Intent     : %d", x$intent_code))
    if (!is.na(x$intent_name) && nzchar(x$intent_name)) {
      cat(sprintf(" (%s)", x$intent_name))
    }
    cat("\n")
  }
  if (!is.na(x$descrip) && nzchar(x$descrip)) {
    cat(sprintf("  Description: %s\n", x$descrip))
  }
  if (!is.na(x$scl_slope) && x$scl_slope != 0 && x$scl_slope != 1) {
    cat(sprintf("  Scaling    : slope=%.4g, intercept=%.4g\n",
                x$scl_slope, x$scl_inter))
  }
  if (!is.na(x$cal_min) && !is.na(x$cal_max) &&
      (x$cal_min != 0 || x$cal_max != 0)) {
    cat(sprintf("  Cal range  : [%.4g, %.4g]\n", x$cal_min, x$cal_max))
  }
  invisible(x)
}
