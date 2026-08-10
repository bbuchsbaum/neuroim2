#' @include common.R
{}
#' @include binary_io.R
{}
#' @include nifti_extensions.R
{}
#' @keywords internal
#' @noRd
.checkDimensions <- function(dimvec) {
  if (!all(dimvec >= 0)) {
    cli::cli_abort("Illegal dimension vector in header: {paste(dimvec, collapse=' x ')}.")
  }
}


## -------------------------------------------------------------------------
## Choosing an output datatype, and the scaling that goes with it
##
## Two things were wrong before. Writing always defaulted to FLOAT, so an
## int16 image doubled in size on a round trip; and asking for an integer type
## went through `writeBin(as.integer(x))`, which truncates toward zero without
## a word of warning -- round-trip error 0.998 on data spanning +/-3.7.
##
## The rules now:
##   * With no `data_type`, keep the source file's datatype when the values are
##     still exactly representable in it (integral, in range, allowing for the
##     source's own slope/intercept). That makes read-write bit-identical for
##     unmodified images. Otherwise fall back to FLOAT, so computed data never
##     silently loses precision to an integer type nobody asked for.
##   * With an explicit integer `data_type`, fit the data to it with
##     scl_slope/scl_inter the way nibabel does, instead of truncating. Data
##     that is already integral and in range is written unscaled, so masks and
##     label volumes stay exact.
## -------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.dtype_limits <- function(data_type) {
  switch(toupper(data_type),
    UBYTE  = c(0, 255),
    BYTE   = c(-128, 127),
    SHORT  = c(-32768, 32767),
    USHORT = c(0, 65535),
    INT    = c(-2147483648, 2147483647),
    UINT   = c(0, 4294967295),
    LONG   = c(-9223372036854775808, 9223372036854774784),
    ULONG  = c(0, 18446744073709549568),
    NULL
  )
}

#' @keywords internal
#' @noRd
.is_integer_dtype <- function(data_type) !is.null(.dtype_limits(data_type))

#' Scaling that maps `vals` into `data_type` without loss where possible.
#' @return list(slope, intercept); the stored value is (x - intercept)/slope
#' @keywords internal
#' @noRd
.nifti_scaling_for <- function(vals, data_type) {
  lim <- .dtype_limits(data_type)
  if (is.null(lim)) return(list(slope = 1, intercept = 0))

  finite <- vals[is.finite(vals)]
  if (length(finite) == 0L) return(list(slope = 1, intercept = 0))

  mn <- min(finite); mx <- max(finite)
  # Exactly representable already? Then do not scale -- masks, label volumes
  # and integer atlases must stay byte-for-byte what they were.
  if (mn >= lim[1] && mx <= lim[2] && all(finite == trunc(finite))) {
    return(list(slope = 1, intercept = 0))
  }
  if (mx == mn) return(list(slope = 1, intercept = mn))

  slope <- (mx - mn) / (lim[2] - lim[1])
  if (!is.finite(slope) || slope == 0) return(list(slope = 1, intercept = 0))
  intercept <- mn - slope * lim[1]
  if (!is.finite(intercept)) return(list(slope = 1, intercept = 0))
  list(slope = slope, intercept = intercept)
}

#' Does the source file's own scaling still represent these values exactly?
#' @keywords internal
#' @noRd
.source_scaling_is_exact <- function(vals, source_header, data_type) {
  lim <- .dtype_limits(data_type)
  if (is.null(lim)) return(FALSE)
  slope <- source_header$scl_slope
  inter <- source_header$scl_intercept
  if (is.null(slope) || !is.finite(slope) || slope == 0) slope <- 1
  if (is.null(inter) || !is.finite(inter)) inter <- 0

  finite <- vals[is.finite(vals)]
  if (length(finite) == 0L) return(TRUE)
  stored <- (finite - inter) / slope
  all(stored >= lim[1]) && all(stored <= lim[2]) && all(stored == round(stored))
}

#' @keywords internal
#' @noRd
.source_scaling <- function(source_header) {
  slope <- source_header$scl_slope
  inter <- source_header$scl_intercept
  if (is.null(slope) || !is.finite(slope) || slope == 0) slope <- 1
  if (is.null(inter) || !is.finite(inter)) inter <- 0
  list(slope = as.numeric(slope), intercept = as.numeric(inter))
}

#' Decide the output datatype and the scaling that goes with it.
#'
#' Keeping the source datatype is only worth doing if the source's *own*
#' slope/intercept still describe the values; re-deriving a fresh slope for the
#' same type would quantise data that was previously exact.
#' @keywords internal
#' @noRd
.nifti_output_plan <- function(vals, source_header, data_type = NULL) {
  src <- source_header$data_storage
  have_src <- !is.null(src) && nzchar(src) && !identical(src, "UNKNOWN")

  if (is.null(data_type)) {
    if (!have_src) {
      data_type <- "FLOAT"
    } else if (!.is_integer_dtype(src)) {
      data_type <- src
    } else if (.source_scaling_is_exact(vals, source_header, src)) {
      return(c(list(data_type = src), .source_scaling(source_header)))
    } else {
      # Computed values no longer fit the file's integer type; promoting keeps
      # them exact rather than quietly quantising them to fit.
      data_type <- "FLOAT"
    }
  } else {
    data_type <- toupper(data_type)
    if (have_src && identical(data_type, src) &&
        .source_scaling_is_exact(vals, source_header, src)) {
      return(c(list(data_type = src), .source_scaling(source_header)))
    }
  }

  c(list(data_type = data_type), .nifti_scaling_for(vals, data_type))
}

#' The voxel values of an image as a plain double vector, without copying when
#' the object already is one.
#'
#' \code{as.vector()} on a dense image duplicates the whole payload -- 58 MB on
#' an ordinary structural, which measured as more than the compiled write it was
#' feeding. An object extending \code{array} already holds exactly the vector we
#' want in \code{.Data}.
#' @keywords internal
#' @noRd
.image_values <- function(x) {
  if (isS4(x) && methods::.hasSlot(x, ".Data")) {
    d <- x@.Data
    if (is.double(d) && !is.object(d)) return(d)
  }
  as.numeric(as.vector(x))
}

#' The voxel values of an image as a plain array, keeping \code{dim}.
#'
#' Same reason as \code{.image_values}, for the callers that need the shape:
#' \code{as.array()} on a dense image goes through S4 dispatch and hands back a
#' fresh copy, which measured as much as the filter kernel it was feeding.
#' @keywords internal
#' @noRd
.image_array <- function(x) {
  if (isS4(x) && methods::.hasSlot(x, ".Data")) {
    d <- x@.Data
    if (is.double(d) && !is.object(d) && !is.null(dim(d))) return(d)
  }
  as.array(x)
}

#' Header bytes plus payload, in one compiled pass.
#' @keywords internal
#' @noRd
.write_nifti_file <- function(hdr, file_name, vals, data_type) {
  raw_header <- .nifti_header_raw(hdr)
  gzipped <- endsWith(file_name, ".gz")
  n <- nifti_write_data_cpp(file_name, raw_header, vals,
                            .getDataCode(data_type),
                            hdr$scl_slope, hdr$scl_intercept,
                            !identical(hdr$endian, .Platform$endian),
                            gzipped)
  invisible(n)
}

#' @keywords internal
#' @noRd
write_nifti_vector <- function(vec, file_name, data_type=NULL, version=NULL) {
  if (length(dim(vec)) != 4) {
    cli::cli_abort("Input vector must be 4-dimensional, not {length(dim(vec))}D.")
  }
	ext <- .make_nifti_volume_labels_extension(volume_labels(vec))
	vals <- .image_values(vec)
	hdr <- as_nifti_header(vec, file_name = file_name, data_type = data_type,
                         extensions = ext, values = vals, version = version)
	.write_nifti_file(hdr, file_name, vals, hdr$data_storage)
}

#' @keywords internal
#' @noRd
write_nifti_hyper_vector <- function(vec, file_name, data_type=NULL, version=NULL) {
  if (!inherits(vec, "NeuroHyperVec")) {
    cli::cli_abort("{.arg vec} must be a {.cls NeuroHyperVec} object.")
  }
  if (length(dim(vec)) != 5) {
    cli::cli_abort("Input hyper-vector must be 5-dimensional, not {length(dim(vec))}D.")
  }

	vals <- as.numeric(dense_array_5d(vec))
	hdr <- as_nifti_header(vec, file_name = file_name, data_type = data_type,
	                       values = vals, version = version)
	# A hyper-vector's 4th and 5th axes are not time, so leave their pixdim at 1.
	hdr$pixdim[5:8] <- 1
	.write_nifti_file(hdr, file_name, vals, hdr$data_storage)
}

#' @keywords internal
#' @noRd
write_nifti_volume <- function(vol, file_name, data_type=NULL, version=NULL) {
  if (length(dim(vol)) != 3) {
    cli::cli_abort("Input volume must be 3-dimensional, not {length(dim(vol))}D.")
  }
	vals <- .image_values(vol)
	hdr <- as_nifti_header(vol, file_name=file_name, data_type=data_type,
	                       values = vals, version = version)
	.write_nifti_file(hdr, file_name, vals, hdr$data_storage)
}


#' Construct a Minimal NIfTI-1 Header from a NeuroVol
#'
#' @description
#' Given a \code{\link[neuroim2]{NeuroVol}} object (or similar), this function
#' builds a basic NIfTI-1 header structure, populating essential fields such as
#' \code{dim}, \code{pixdim}, \code{datatype}, the affine transform, and the
#' quaternion parameters.
#'
#' @details
#' This is a convenience function that calls \code{\link{createNIfTIHeader}}
#' first, then updates the fields (dimensions, \code{pixdim}, orientation, etc.)
#' based on the \code{vol} argument. The voxel offset is set to 352 bytes (or
#' larger if extensions are provided), and the quaternion is derived from the
#' transform matrix via \code{\link{matrixToQuatern}}.
#'
#' Note: This function primarily sets up a minimal header suitable for writing
#' standard single-file NIfTI-1. If you need a more comprehensive or advanced
#' usage, consider manually editing the returned list.
#'
#' @param vol A \code{\link[neuroim2]{NeuroVol}} (or 3D array-like) specifying
#'   dimensions, spacing, and affine transform.
#' @param file_name A character string for the file name (used within the header
#'   but not necessarily to write data).
#' @param oneFile Logical; if \code{TRUE}, sets the NIfTI magic to \code{"n+1"},
#'   implying a single-file format (\code{.nii}). If \code{FALSE}, uses
#'   \code{"ni1"} (header+image).
#' @param data_type Character specifying the data representation, e.g. \code{"FLOAT"},
#'   \code{"DOUBLE"}, \code{"SHORT"}. \code{NULL} (the default) keeps the source
#'   file's datatype when \code{vol} came from one and its values are still
#'   exactly representable in it, and uses \code{"FLOAT"} otherwise.
#' @param extensions Optional \code{\link{NiftiExtensionList-class}} object or list
#'   of \code{\link{NiftiExtension-class}} objects to include in the header.
#' @param values Optional numeric vector of the voxel values that will be
#'   written. Used to pick a datatype and to derive \code{scl_slope} /
#'   \code{scl_inter} for integer output; taken from \code{vol} when omitted.
#' @param version Either \code{1} or \code{2}, selecting the NIfTI version to
#'   describe. \code{NULL} (the default) uses NIfTI-1 unless a dimension exceeds
#'   what its 16-bit \code{dim} field can hold, in which case NIfTI-2 is used
#'   because NIfTI-1 physically cannot represent the image.
#'
#' @return A \code{list} representing the NIfTI header fields, containing
#'   elements like \code{dimensions}, \code{pixdim}, \code{datatype}, \code{qform},
#'   \code{quaternion}, \code{qfac}, \code{extensions}, etc. This can be passed to
#'   other functions that write or manipulate the header.
#'
#' @details
#' When \code{vol} was read from a NIfTI file, the fields that describe the
#' acquisition rather than the array -- repetition time (\code{pixdim[5]}),
#' \code{xyzt_units}, \code{qform_code} and \code{sform_code}, \code{descrip},
#' \code{aux_file}, the intent fields, \code{cal_min}/\code{cal_max}, the slice
#' timing fields and \code{toffset} -- are carried over from that file. Only the
#' fields the object itself determines (geometry, dimensions, datatype, offsets)
#' are recomputed. Before this behaviour existed, a read-write round trip reset
#' all of them, which silently dropped the TR and relabelled an MNI-space image
#' as scanner-space.
#'
#' @seealso
#' \code{\link{createNIfTIHeader}} for the base constructor of an empty NIfTI header.
#' \code{\link{NiftiExtension}} for creating extensions.
#'
#' @export
as_nifti_header <- function(vol, file_name, oneFile=TRUE, data_type=NULL,
                           extensions = NULL, values = NULL, version = NULL) {
		src <- .source_nifti_header(vol)

		# Start from the source header when there is one, so acquisition
		# metadata survives; createNIfTIHeader() supplies neutral defaults
		# otherwise.
		hd <- createNIfTIHeader(oneFile=oneFile, file_name=file_name)
		if (length(src)) {
		  carry <- c("diminfo", "intent1", "intent2", "intent3", "intent_code",
		             "intent_name", "slice_start", "slice_end", "slice_code",
		             "slice_duration", "toffset", "xyzt_units",
		             "cal_max", "cal_min", "description", "auxfile",
		             "qform_code", "sform_code")
		  for (nm in carry) {
		    if (!is.null(src[[nm]])) hd[[nm]] <- src[[nm]]
		  }
		}

		if (is.null(values)) values <- .image_values(vol)

		plan <- .nifti_output_plan(values, src, data_type)
		hd$datatype <- .getDataCode(plan$data_type)
		hd$data_storage <- .getDataStorage(hd$datatype)
		hd$bitpix <- .getDataSize(hd$data_storage) * 8

		hd$file_name <- file_name
		hd$endian <- .Platform$endian
		hd$dimensions <- c(length(dim(vol)), dim(vol))
		N <- 8 - length(hd$dimensions)
		hd$dimensions <- c(hd$dimensions,  rep(1, N))
		hd$num_dimensions <- length(dim(vol))

		# pixdim[2:4] is geometry and comes from the object; pixdim[5:8] is
		# timing and beyond, which only the source file knows.
		hd$pixdim <- c(0, spacing(vol), rep(0, 4))
		if (length(src$pixdim) == 8L) {
		  keep <- src$pixdim[5:8]
		  keep[!is.finite(keep)] <- 0
		  hd$pixdim[5:8] <- keep
		}

		hd$scl_slope <- plan$slope
		hd$scl_intercept <- plan$intercept

		tmat <- trans(vol)

		# Derive qoffset from the transform matrix translation column
		# to guarantee consistency between qoffset and sform
		hd$qoffset <- tmat[1:3, 4]

		hd$qform <- tmat
		hd$sform <- tmat

		quat1 <- matrixToQuatern(tmat)
		hd$quaternion <- quat1$quaternion
		hd$qfac <- quat1$qfac
		hd$pixdim[1] <- hd$qfac

		# A code of 0 means "no transform here", which would discard the affine
		# we just wrote. Only meaningful when the source really had none.
		if (!is.numeric(hd$qform_code) || hd$qform_code < 1) hd$qform_code <- 1
		if (!is.numeric(hd$sform_code) || hd$sform_code < 1) hd$sform_code <- 1

		# Handle extensions
		if (!is.null(extensions)) {
		  if (is.list(extensions) && !is(extensions, "NiftiExtensionList")) {
		    # Convert list of NiftiExtension objects to NiftiExtensionList
		    extensions <- new("NiftiExtensionList", extensions)
		  }
		  hd$extensions <- extensions
		} else {
		  hd$extensions <- new("NiftiExtensionList")
		}

		hd$version <- .nifti_output_version(hd$dimensions, version)
		hd$magic <- if (hd$version == 2L) {
		  if (oneFile) "n+2" else "ni2"
		} else {
		  if (oneFile) "n+1" else "ni1"
		}
		hd$onefile <- oneFile

		# vox_offset: fixed header size plus extension bytes.
		# total_extension_size includes the 4-byte extender field.
		hd$vox_offset <- .nifti_header_size(hd$version) +
		  total_extension_size(hd$extensions)

		hd

}

#' The raw NIfTI header an image was read from, or an empty list.
#' @keywords internal
#' @noRd
.source_nifti_header <- function(x) {
  hdr <- tryCatch(
    if (isS4(x) && methods::.hasSlot(x, "header")) methods::slot(x, "header") else NULL,
    error = function(e) NULL
  )
  if (is.list(hdr) && length(hdr)) hdr else list()
}

#' Record the file header an object was read from.
#' @keywords internal
#' @noRd
.attach_header <- function(obj, meta) {
  hdr <- tryCatch(
    if (isS4(meta) && methods::.hasSlot(meta, "header")) methods::slot(meta, "header") else NULL,
    error = function(e) NULL
  )
  if (is.list(hdr) && length(hdr) && isS4(obj) && methods::.hasSlot(obj, "header")) {
    obj@header <- hdr
  }
  obj
}

#' @keywords internal
#' @noRd
.nifti_header_size <- function(version) if (identical(as.integer(version), 2L)) 540L else 348L

#' @keywords internal
#' @noRd
.nifti_output_version <- function(dimensions, version = NULL) {
  if (!is.null(version)) {
    v <- as.integer(version)
    if (!v %in% c(1L, 2L)) {
      cli::cli_abort("{.arg version} must be 1 or 2, not {.val {version}}.")
    }
    if (v == 1L && any(dimensions > 32767)) {
      cli::cli_abort(c(
        "NIfTI-1 cannot store a dimension larger than 32767.",
        "x" = "This image is {paste(dimensions[-1][seq_len(dimensions[1])],
                                    collapse = ' x ')}.",
        "i" = "Use {.code version = 2}."
      ))
    }
    return(v)
  }
  if (any(dimensions > 32767)) 2L else 1L
}


#' Create an Empty NIfTI-1 Header List
#'
#' @description
#' Initializes a list of fields following the NIfTI-1 specification with default
#' or placeholder values. Users typically call this internally via
#' \code{\link{as_nifti_header}} rather than using directly.
#'
#' @details
#' This function sets up the skeleton of a NIfTI-1 header, including fields for
#' \code{diminfo}, \code{pixdim}, \code{qform_code}, \code{magic}, etc. Most fields
#' are initialized to zero, empty characters, or standard placeholders. The
#' \code{oneFile} argument controls whether \code{"n+1"} or \code{"ni1"} is used
#' for the \code{magic} field.
#'
#' @param oneFile Logical; if \code{TRUE}, \code{magic} is set to \code{"n+1"}
#'   indicating a single-file (.nii) approach. Otherwise set to \code{"ni1"}.
#' @param file_name Optional character string to store in the header, usually
#'   referencing the intended output file name.
#'
#' @return A named \code{list} containing approximately 30 fields that comprise
#'   the NIfTI-1 header structure. Many of these are placeholders until filled by
#'   downstream usage.
#'
#' @seealso
#' \code{\link{as_nifti_header}} for populating the returned header with actual
#' data from a \code{NeuroVol}.
#'
#' @export
createNIfTIHeader <- function(oneFile=TRUE, file_name=NULL) {
	header <- list()
	header$file_type <- "NIfTI"
	header$encoding <- "binary"
	header$version <- "1"


	header$file_name <- file_name
	header$endian <- .Platform$endian

	header$diminfo <- 0
	header$dimensions <- NULL
	header$num_dimensions <- NULL

	header$intent1 <-  0
	header$intent2 <-  0
	header$intent3 <-  0


	header$intent_code <-  0
	header$datatype <- NULL
	header$data_storage <- NULL
	header$bitpix <- NULL
	header$slice_start <- 0
	header$pixdim <-  NULL

	header$qfac <- -1


	header$vox_offset <- 0
	header$scl_slope <- 0
	header$scl_intercept <- 0
	header$slice_end <- 0
	header$slice_code <- 0


	header$xyzt_units <- 2
	header$cal_max <- 0
	header$cal_min <- 0

	header$slice_duration <- 0
	header$toffset <- 0

	header$glmax <- 0
	header$glmin <- 0

	header$description <- integer(80)
	header$auxfile <- integer(24)

	header$qform_code <- 1
	header$sform_code <- 1
	header$quaternion <- NULL

	header$qoffset <- NULL
	header$qform <- NULL

	header$sform <- NULL
	header$intent_name <- integer(16)
	header$magic <- "n+1"

	header$onefile <- oneFile
	if (oneFile) {
		header$magic <- "n+1"
	} else {
		header$magic <- "ni1"
	}

	header$version <- 1
	header

}


#' Read NIfTI Header
#'
#' @description
#' Reads a NIfTI-1 header from file, including any extensions present.
#'
#' @param fname Character string specifying the path to the NIfTI file
#'   (.nii, .nii.gz, or .hdr).
#' @param read_extensions Logical; if TRUE (default), reads any extensions
#'   present in the file. Set to FALSE to skip extension reading.
#'
#' @return A list containing header fields including:
#'   \itemize{
#'     \item Standard NIfTI header fields (dimensions, pixdim, transforms, etc.)
#'     \item \code{extensions}: A \code{\link{NiftiExtensionList-class}} object
#'       containing any extensions found (empty if none present or read_extensions=FALSE)
#'   }
#'
#' The NIfTI spec's METHOD 1 affine: voxel size only, image centred on the
#' origin, with the first axis flipped. This is what a file that declares
#' \code{qform_code == 0} and \code{sform_code == 0} is asking for -- both
#' transforms are marked "unknown", so the quaternion in the header is not a
#' statement about where the image is.
#' @keywords internal
#' @noRd
.method1_affine <- function(shape, zooms) {
  shape <- as.numeric(shape); zooms <- as.numeric(zooms)
  n <- min(3L, length(shape), length(zooms))
  z <- rep(1, 3); z[seq_len(n)] <- zooms[seq_len(n)]
  s <- rep(1, 3); s[seq_len(n)] <- shape[seq_len(n)]
  z[!is.finite(z) | z == 0] <- 1
  z[1] <- -abs(z[1])                    # NIfTI is RAS+, arrays come in LAS
  z[2:3] <- abs(z[2:3])
  m <- diag(c(z, 1))
  m[1:3, 4] <- -(s - 1) / 2 * z
  m
}

#' Read 8-byte little/big-endian integers as doubles.
#'
#' R has no 64-bit integer type, and every NIfTI-2 field that uses one (dims,
#' vox_offset, slice indices) is far below 2^53, so a double holds them exactly.
#' @keywords internal
#' @noRd
.read_int64 <- function(conn, n, endian) {
  raw_bytes <- readBin(conn, "raw", n = n * 8L)
  if (length(raw_bytes) != n * 8L) {
    stop("Unexpected end of file while reading the NIfTI-2 header")
  }
  vapply(seq_len(n), function(i) {
    b <- as.integer(raw_bytes[((i - 1L) * 8L + 1L):(i * 8L)])
    if (endian == "big") b <- rev(b)
    val <- sum(b * 256^(0:7))
    if (b[8] >= 128) val <- val - 2^64   # two's complement
    val
  }, numeric(1))
}

#' @keywords internal
#' @noRd
.open_nifti_conn <- function(fname) {
  if (.isExtension(fname, ".nii.gz") || .isExtension(fname, ".hdr.gz") ||
      endsWith(tolower(fname), ".gz")) {
    list(conn = gzfile(fname, open = "rb"), encoding = "gzip")
  } else if (.isExtension(fname, ".nii") || .isExtension(fname, ".hdr")) {
    list(conn = file(fname, open = "rb"), encoding = "binary")
  } else {
    stop(paste("illegal NIFTI header name", fname))
  }
}

#' Byte order and version, from the leading sizeof_hdr field.
#' @keywords internal
#' @noRd
.nifti_endian_and_version <- function(conn, fname) {
  first <- readBin(conn, "raw", n = 4L)
  if (length(first) < 4L) {
    cli::cli_abort(c(
      "{.path {fname}} is not a NIfTI or ANALYZE file.",
      "x" = "It is only {length(first)} byte{?s} long; a header needs at least 348."
    ))
  }
  for (endian in c("little", "big")) {
    size <- readBin(first, integer(), n = 1L, size = 4L, endian = endian)
    if (identical(size, 348L)) return(list(endian = endian, version = 1L))
    if (identical(size, 540L)) return(list(endian = endian, version = 2L))
  }
  le <- readBin(first, integer(), n = 1L, size = 4L, endian = "little")
  cli::cli_abort(c(
    "{.path {fname}} does not start with a NIfTI or ANALYZE header.",
    "x" = "The leading {.field sizeof_hdr} field is {.val {le}}; it must be
           348 (NIfTI-1/ANALYZE) or 540 (NIfTI-2), in either byte order.",
    "i" = "The file may be compressed with something other than gzip, be
           truncated, or not be an image at all."
  ))
}

#' @keywords internal
#' @noRd
read_nifti_header <- function(fname, read_extensions = TRUE) {
	header <- list()
	header$file_type <- "NIfTI"

	opened <- .open_nifti_conn(fname)
	conn <- opened$conn
	header$encoding <- opened$encoding

	# Ensure connection is closed on exit
	on.exit(close(conn), add = TRUE)

	ev <- .nifti_endian_and_version(conn, fname)
	endian <- ev$endian

	header$file_name <- fname
	header$endian <- endian

	header <- if (ev$version == 2L) {
	  .read_nifti2_fields(conn, header, endian)
	} else {
	  .read_nifti1_fields(conn, header, endian)
	}

	header$data_storage <- .getDataStorage(header$datatype)

	# ANALYZE 7.5 has no magic and no transform fields; everything from
	# qform_code onwards is a different struct there, so trust only the
	# geometry and derive the affine the way the reference implementations do.
	header$is_analyze <- !(header$magic %in% c("n+1", "ni1", "n+2", "ni2"))

	header$qfac <- header$pixdim[1]
	if (!is.finite(header$qfac) || header$qfac == 0) header$qfac <- 1

	shape <- header$dimensions[-1][seq_len(max(1L, header$num_dimensions))]
	if (header$is_analyze) {
	  header$qform_code <- 0L
	  header$sform_code <- 0L
	  header$qform <- .method1_affine(shape, header$pixdim[2:4])
	  header$sform <- header$qform
	} else {
	  header$qform <- quaternToMatrix(header$quaternion, header$qoffset,
	                                  header$pixdim[2:4], header$qfac)
	  header$sform <- rbind(matrix(header$srow, 3, 4, byrow = TRUE), c(0, 0, 0, 1))
	  if (header$qform_code <= 0 && header$sform_code <= 0) {
	    # Both transforms are declared unknown: METHOD 1.
	    header$qform <- .method1_affine(shape, header$pixdim[2:4])
	    header$sform <- header$qform
	  }
	}

	header$onefile <- substr(header$magic, 2, 2) == "+"
	header$version <- ev$version

	# Extensions live between the fixed header and vox_offset, and only in the
	# single-file layout.
	hdr_size <- .nifti_header_size(ev$version)
	if (read_extensions && header$onefile && header$vox_offset > hdr_size + 4) {
		header$extensions <- read_nifti_extensions(conn, header$vox_offset, endian)
	} else {
		header$extensions <- new("NiftiExtensionList")
	}

	header

}

#' @keywords internal
#' @noRd
.read_nifti1_fields <- function(conn, header, endian) {
	readBin(conn, what=integer(), n=10+18+4+2+1, size=1)

	header$diminfo <- readBin(conn, what=integer(), n=1, size=1)
	header$dimensions <- readBin(conn, integer(), n=8, size=2, endian=endian)
	header <- .normalize_nifti_dims(header)

	header$intent1 <-  readBin(conn, double(), n=1, size=4, endian=endian)
	header$intent2 <-  readBin(conn, double(), n=1, size=4, endian=endian)
	header$intent3 <-  readBin(conn, double(), n=1, size=4, endian=endian)

	header$intent_code <-  readBin(conn, integer(), n=1, size=2, endian=endian)
	header$datatype <- readBin(conn, integer(), n=1, size=2, endian=endian)
	header$bitpix <- readBin(conn, integer(), n=1, size=2, endian=endian)
	header$slice_start <- readBin(conn, integer(), n=1, size=2, endian=endian)
	header$pixdim <-  readBin(conn, double(), n=8, size=4, endian=endian)

	header$vox_offset <- readBin(conn, double(), n=1, size=4, endian=endian)
	header$scl_slope <- readBin(conn, double(), n=1, size=4, endian=endian)
	header$scl_intercept <- readBin(conn, double(), n=1, size=4, endian=endian)
	header$slice_end <- readBin(conn, integer(), n=1, size=2, endian=endian)
	header$slice_code <-  readBin(conn, integer(), n=1, size=1, endian=endian)

	header$xyzt_units <- readBin(conn, integer(), n=1, size=1, endian=endian)
	header$cal_max <- readBin(conn, double(), n=1, size=4, endian=endian)
	header$cal_min <- readBin(conn, double(), n=1, size=4, endian=endian)

	header$slice_duration <- readBin(conn, double(), n=1, size=4, endian=endian)
	header$toffset <- readBin(conn, double(), n=1, size=4, endian=endian)

	header$glmax <- readBin(conn, integer(), n=1, size=4, endian=endian) # unused glmax, glmin
	header$glmin <- readBin(conn, integer(), n=1, size=4, endian=endian) # unused glmax, glmin

	header$description <- readBin(conn, integer(), n=80, size=1, signed=FALSE)
	header$auxfile <- readBin(conn, integer(), n=24, size=1, signed=FALSE)

	header$qform_code <- readBin(conn, integer(), n=1, size=2, endian=endian)
	header$sform_code <- readBin(conn, integer(), n=1, size=2, endian=endian)
	header$quaternion <- readBin(conn, double(), n=3, size=4, endian=endian)

	header$qoffset <- readBin(conn, double(), n=3, size=4, endian=endian)
	header$srow  <- readBin(conn, double(), n=12, size=4, endian=endian)
	# Fixed-width byte field, not a C string: reading it as `character()` used
	# to run past the field whenever the 16 bytes held no null terminator.
	header$intent_name <- readBin(conn, integer(), n=16, size=1, signed=FALSE)
	header$magic <- .bytes_to_string(readBin(conn, integer(), n=4, size=1, signed=FALSE))

	header
}

#' @keywords internal
#' @noRd
.read_nifti2_fields <- function(conn, header, endian) {
	header$magic <- .bytes_to_string(readBin(conn, integer(), n=8, size=1, signed=FALSE))
	header$datatype <- readBin(conn, integer(), n=1, size=2, endian=endian)
	header$bitpix <- readBin(conn, integer(), n=1, size=2, endian=endian)
	header$dimensions <- .read_int64(conn, 8L, endian)
	header <- .normalize_nifti_dims(header)

	header$intent1 <- readBin(conn, double(), n=1, size=8, endian=endian)
	header$intent2 <- readBin(conn, double(), n=1, size=8, endian=endian)
	header$intent3 <- readBin(conn, double(), n=1, size=8, endian=endian)
	header$pixdim <- readBin(conn, double(), n=8, size=8, endian=endian)
	header$vox_offset <- .read_int64(conn, 1L, endian)
	header$scl_slope <- readBin(conn, double(), n=1, size=8, endian=endian)
	header$scl_intercept <- readBin(conn, double(), n=1, size=8, endian=endian)
	header$cal_max <- readBin(conn, double(), n=1, size=8, endian=endian)
	header$cal_min <- readBin(conn, double(), n=1, size=8, endian=endian)
	header$slice_duration <- readBin(conn, double(), n=1, size=8, endian=endian)
	header$toffset <- readBin(conn, double(), n=1, size=8, endian=endian)
	header$slice_start <- .read_int64(conn, 1L, endian)
	header$slice_end <- .read_int64(conn, 1L, endian)
	header$description <- readBin(conn, integer(), n=80, size=1, signed=FALSE)
	header$auxfile <- readBin(conn, integer(), n=24, size=1, signed=FALSE)
	header$qform_code <- readBin(conn, integer(), n=1, size=4, endian=endian)
	header$sform_code <- readBin(conn, integer(), n=1, size=4, endian=endian)
	header$quaternion <- readBin(conn, double(), n=3, size=8, endian=endian)
	header$qoffset <- readBin(conn, double(), n=3, size=8, endian=endian)
	header$srow <- readBin(conn, double(), n=12, size=8, endian=endian)
	header$slice_code <- readBin(conn, integer(), n=1, size=4, endian=endian)
	header$xyzt_units <- readBin(conn, integer(), n=1, size=4, endian=endian)
	header$intent_code <- readBin(conn, integer(), n=1, size=4, endian=endian)
	header$intent_name <- readBin(conn, integer(), n=16, size=1, signed=FALSE)
	header$diminfo <- readBin(conn, integer(), n=1, size=1, signed=FALSE)
	readBin(conn, integer(), n=15, size=1)   # unused_str

	header$glmax <- 0
	header$glmin <- 0
	header
}

#' Normalise the `dim` array without inventing or discarding axes.
#'
#' The old rule took dimensions up to the first entry equal to 1, which threw
#' away a trailing singleton (a 1-volume run became a 3-D volume) and, worse,
#' truncated a genuine single-slice image: 64 x 64 x 1 came back as 64 x 64.
#' The standard says dim[0] is the count; use it.
#' @keywords internal
#' @noRd
.normalize_nifti_dims <- function(header) {
  d <- header$dimensions
  nd <- d[1]
  if (!is.finite(nd) || nd < 1 || nd > 7) {
    cli::cli_abort(c(
      "Invalid NIfTI header: {.field dim[0]} is {.val {nd}}.",
      "i" = "It must be between 1 and 7."
    ))
  }
  nd <- as.integer(nd)
  # Only the axes dim[0] declares are meaningful; the rest are padding and are
  # normalised to 1 so downstream arithmetic (prod(), offsets) stays safe.
  shape <- d[2:(nd + 1L)]
  shape[!is.finite(shape) | shape == 0] <- 1
  .checkDimensions(shape)
  header$dimensions <- c(nd, shape, rep(1, 7L - nd))
  header$num_dimensions <- nd
  header
}

#' @keywords internal
#' @noRd
.bytes_to_string <- function(b) {
  b <- as.integer(b)
  stop_at <- which(b == 0)[1]
  if (!is.na(stop_at)) b <- b[seq_len(stop_at - 1L)]
  if (!length(b)) return("")
  rawToChar(as.raw(b))
}


#' Write NIfTI Header
#'
#' @description
#' Writes a NIfTI-1 header to a connection, including any extensions.
#'
#' @param niftiInfo A list containing NIfTI header fields.
#' @param conn An open binary connection to write to.
#' @param close Logical; if TRUE (default), close the connection after writing.
#'
#' @return The connection (invisibly).
#'
#' @details
#' The function writes the 348-byte standard header, then extensions (if present),
#' then pads to vox_offset. The vox_offset in niftiInfo should already account
#' for extension size.
#'
#' @keywords internal
#' @noRd
write_nifti_header <- function(niftiInfo, conn, close=TRUE) {
	writeBin(.nifti_header_raw(niftiInfo), conn)
	if (close) {
		close(conn)
	}
	invisible(conn)
}

#' Coerce a header text field to exactly `n` zero-padded bytes.
#'
#' Accepts either a character string or the byte vector a reader produced, and
#' truncates rather than overrunning the fixed-width field.
#' @keywords internal
#' @noRd
.string_to_bytes <- function(x, n) {
  b <- if (is.character(x)) {
    s <- paste(x[nzchar(x)], collapse = "")
    if (nzchar(s)) as.integer(charToRaw(s)) else integer(0)
  } else if (is.raw(x)) {
    as.integer(x)
  } else if (is.numeric(x)) {
    v <- as.integer(x)
    v[is.na(v)] <- 0L
    v
  } else {
    integer(0)
  }
  b <- b[seq_len(min(length(b), n))]
  c(b, rep(0L, n - length(b)))
}

#' @keywords internal
#' @noRd
.hdr_num <- function(x, default = 0) {
  if (is.null(x) || length(x) == 0L) return(default)
  x <- as.numeric(x)[1]
  if (!is.finite(x)) default else x
}

#' Serialize a header list to the bytes that go at the front of the file.
#'
#' Producing bytes rather than writing to a connection is what lets the header
#' and the payload be handed to the compiled writer in one call, and it is also
#' the only way to build a gzipped file as a single stream.
#'
#' @return a raw vector of length \code{vox_offset}
#' @keywords internal
#' @noRd
.nifti_header_raw <- function(niftiInfo) {
  version <- if (is.null(niftiInfo$version)) 1L else as.integer(niftiInfo$version)
  if (!version %in% c(1L, 2L)) version <- 1L
  endian <- if (is.null(niftiInfo$endian)) .Platform$endian else niftiInfo$endian

  stopifnot(length(niftiInfo$dimensions) == 8)
  stopifnot(length(niftiInfo$pixdim) == 8)

  con <- rawConnection(raw(0), "wb")
  on.exit(close(con), add = TRUE)

  sform <- niftiInfo$sform
  srow <- if (is.matrix(sform)) as.double(t(sform[1:3, , drop = FALSE])) else rep(0, 12)

  if (version == 1L) {
    writeBin(348L, con, 4, endian)
    writeBin(integer(34), con, 1, endian)          # data_type[10], db_name[18], extents, session_error
    writeBin(charToRaw("r"), con)                  # regular
    writeBin(as.integer(.hdr_num(niftiInfo$diminfo)), con, size = 1, endian)
    writeBin(as.integer(niftiInfo$dimensions), con, 2, endian)
    writeBin(as.double(.hdr_num(niftiInfo$intent1)), con, 4, endian)
    writeBin(as.double(.hdr_num(niftiInfo$intent2)), con, 4, endian)
    writeBin(as.double(.hdr_num(niftiInfo$intent3)), con, 4, endian)
    writeBin(as.integer(.hdr_num(niftiInfo$intent_code)), con, 2, endian)
    writeBin(as.integer(niftiInfo$datatype), con, 2, endian)
    writeBin(as.integer(niftiInfo$bitpix), con, 2, endian)
    writeBin(as.integer(.hdr_num(niftiInfo$slice_start)), con, 2, endian)
    writeBin(as.double(niftiInfo$pixdim), con, 4, endian)
    writeBin(as.double(niftiInfo$vox_offset), con, 4, endian)
    writeBin(as.double(.hdr_num(niftiInfo$scl_slope, 1)), con, 4, endian)
    writeBin(as.double(.hdr_num(niftiInfo$scl_intercept)), con, 4, endian)
    writeBin(as.integer(.hdr_num(niftiInfo$slice_end)), con, 2, endian)
    writeBin(as.integer(.hdr_num(niftiInfo$slice_code)), con, 1, endian)
    writeBin(as.integer(.hdr_num(niftiInfo$xyzt_units)), con, 1, endian)
    writeBin(as.double(.hdr_num(niftiInfo$cal_max)), con, 4, endian)
    writeBin(as.double(.hdr_num(niftiInfo$cal_min)), con, 4, endian)
    writeBin(as.double(.hdr_num(niftiInfo$slice_duration)), con, 4, endian)
    writeBin(as.double(.hdr_num(niftiInfo$toffset)), con, 4, endian)
    writeBin(as.integer(.hdr_num(niftiInfo$glmax)), con, 4, endian)
    writeBin(as.integer(.hdr_num(niftiInfo$glmin)), con, 4, endian)
    writeBin(.string_to_bytes(niftiInfo$description, 80L), con, 1, endian)
    writeBin(.string_to_bytes(niftiInfo$auxfile, 24L), con, 1, endian)
    writeBin(as.integer(.hdr_num(niftiInfo$qform_code)), con, 2, endian)
    writeBin(as.integer(.hdr_num(niftiInfo$sform_code)), con, 2, endian)
    writeBin(as.double(niftiInfo$quaternion), con, 4, endian)
    writeBin(as.double(niftiInfo$qoffset), con, 4, endian)
    writeBin(srow, con, 4, endian)
    writeBin(.string_to_bytes(niftiInfo$intent_name, 16L), con, 1, endian)
    # magic is 4 bytes: 3 characters and a terminating nul
    writeBin(.string_to_bytes(substr(niftiInfo$magic, 1L, 3L), 4L), con, 1, endian)
  } else {
    writeBin(540L, con, 4, endian)
    # magic is 8 bytes in NIfTI-2 and carries a newline/EOF guard that detects
    # a file mangled by an ASCII-mode transfer.
    magic <- c(.string_to_bytes(niftiInfo$magic, 4L), 0x0D, 0x0A, 0x1A, 0x0A)
    writeBin(as.integer(magic), con, 1, endian)
    writeBin(as.integer(niftiInfo$datatype), con, 2, endian)
    writeBin(as.integer(niftiInfo$bitpix), con, 2, endian)
    .write_int64(con, niftiInfo$dimensions, endian)
    writeBin(as.double(.hdr_num(niftiInfo$intent1)), con, 8, endian)
    writeBin(as.double(.hdr_num(niftiInfo$intent2)), con, 8, endian)
    writeBin(as.double(.hdr_num(niftiInfo$intent3)), con, 8, endian)
    writeBin(as.double(niftiInfo$pixdim), con, 8, endian)
    .write_int64(con, niftiInfo$vox_offset, endian)
    writeBin(as.double(.hdr_num(niftiInfo$scl_slope, 1)), con, 8, endian)
    writeBin(as.double(.hdr_num(niftiInfo$scl_intercept)), con, 8, endian)
    writeBin(as.double(.hdr_num(niftiInfo$cal_max)), con, 8, endian)
    writeBin(as.double(.hdr_num(niftiInfo$cal_min)), con, 8, endian)
    writeBin(as.double(.hdr_num(niftiInfo$slice_duration)), con, 8, endian)
    writeBin(as.double(.hdr_num(niftiInfo$toffset)), con, 8, endian)
    .write_int64(con, .hdr_num(niftiInfo$slice_start), endian)
    .write_int64(con, .hdr_num(niftiInfo$slice_end), endian)
    writeBin(.string_to_bytes(niftiInfo$description, 80L), con, 1, endian)
    writeBin(.string_to_bytes(niftiInfo$auxfile, 24L), con, 1, endian)
    writeBin(as.integer(.hdr_num(niftiInfo$qform_code)), con, 4, endian)
    writeBin(as.integer(.hdr_num(niftiInfo$sform_code)), con, 4, endian)
    writeBin(as.double(niftiInfo$quaternion), con, 8, endian)
    writeBin(as.double(niftiInfo$qoffset), con, 8, endian)
    writeBin(srow, con, 8, endian)
    writeBin(as.integer(.hdr_num(niftiInfo$slice_code)), con, 4, endian)
    writeBin(as.integer(.hdr_num(niftiInfo$xyzt_units)), con, 4, endian)
    writeBin(as.integer(.hdr_num(niftiInfo$intent_code)), con, 4, endian)
    writeBin(.string_to_bytes(niftiInfo$intent_name, 16L), con, 1, endian)
    writeBin(as.integer(.hdr_num(niftiInfo$diminfo)), con, 1, endian)
    writeBin(integer(15), con, 1, endian)          # unused_str
  }

  write_nifti_extensions(con, niftiInfo$extensions, endian)

  out <- rawConnectionValue(con)
  want <- niftiInfo$vox_offset
  if (length(out) > want) {
    cli::cli_abort("Header serialisation produced {length(out)} bytes but
                    {.field vox_offset} is {want}.")
  }
  c(out, raw(want - length(out)))
}

#' Write 64-bit integers as raw bytes; see .read_int64 for why doubles suffice.
#' @keywords internal
#' @noRd
.write_int64 <- function(con, values, endian) {
  for (v in as.numeric(values)) {
    if (!is.finite(v)) v <- 0
    neg <- v < 0
    if (neg) v <- v + 2^64
    b <- integer(8)
    for (i in 1:8) {
      b[i] <- as.integer(v %% 256)
      v <- floor(v / 256)
    }
    if (endian == "big") b <- rev(b)
    writeBin(b, con, 1, endian)
  }
  invisible(NULL)
}
