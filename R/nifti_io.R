#' @include common.R
{}
#' @include binary_io.R
{}

#' @keywords internal
#' @noRd
.checkDimensions <- function(dimvec) {
	if (any(dimvec < 0)) {
		stop(paste("nifti(checkDimensons): illegal dimension vector in header: ", dimvec))
	}
}


#' @keywords internal
#' @noRd
write_nifti_vector <- function(vec, file_name, data_type="FLOAT") {
	assertthat::assert_that(length(dim(vec)) == 4)
	hdr <- as_nifti_header(vec, file_name=file_name, data_type=data_type)

	conn <- if (substr(file_name, nchar(file_name)-2, nchar(file_name)) == ".gz") {
				gzfile(file_name, open="wb")
			} else {
				file(file_name, open="wb")
			}

	write_nifti_header(hdr, conn, close=FALSE)
	writer <- BinaryWriter(conn, hdr$vox_offset, data_type, hdr$bitpix/8, .Platform$endian)

	NVOLS <- dim(vec)[4]

	for (i in 1:NVOLS) {
		write_elements(writer, as.numeric(vec[[i]]))
	}
	close(writer)
}

#' @keywords internal
#' @noRd
write_nifti_volume <- function(vol, file_name, data_type="FLOAT") {
	stopifnot(length(dim(vol)) == 3)
	hdr <- as_nifti_header(vol, file_name=file_name, data_type=data_type)

	conn <- if (substr(file_name, nchar(file_name)-2, nchar(file_name)) == ".gz") {
		gzfile(file_name, open="wb")
	} else {
		file(file_name, open="wb")
	}

	write_nifti_header(hdr, conn, close=FALSE)
	writer <- BinaryWriter(conn, hdr$vox_offset, data_type, hdr$bitpix/8, .Platform$endian)
	write_elements(writer, as.numeric(vol))
	close(writer)
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
#' based on the \code{vol} argument. The voxel offset is set to 352 bytes, and
#' the quaternion is derived from the transform matrix via
#' \code{\link{matrixToQuatern}}.
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
#'   \code{"DOUBLE"}. The internal code picks an integer NIfTI code.
#'
#' @return A \code{list} representing the NIfTI-1 header fields, containing
#'   elements like \code{dimensions}, \code{pixdim}, \code{datatype}, \code{qform},
#'   \code{quaternion}, \code{qfac}, etc. This can be passed to other
#'   functions that write or manipulate the header.
#'
#' @seealso
#' \code{\link{createNIfTIHeader}} for the base constructor of an empty NIfTI header.
#'
#' @export
as_nifti_header <- function(vol, file_name, oneFile=TRUE, data_type="FLOAT") {
		hd <- createNIfTIHeader(oneFile=oneFile, file_name=file_name)
		hd$file_name <- file_name
		hd$endian <- .Platform$endian
		hd$vox_offset <- 352
		hd$datatype <- .getDataCode(data_type)
		hd$data_storage <- .getDataStorage(hd$datatype)
		hd$bitpix <- .getDataSize(data_type) * 8
		hd$dimensions <- c(length(dim(vol)), dim(vol))
		N <- 8 - length(hd$dimensions)
		hd$dimensions <- c(hd$dimensions,  rep(1, N))
		hd$num_dimensions <- length(dim(vol))

		### only encodes pixdim for three dimensions
		hd$pixdim <- c(0, spacing(vol), rep(0,4))

		hd$qoffset <- origin(space(vol))
		hd$scl_intercept <- 0
		hd$scl_slope <- 1

		tmat <- trans(vol)

		hd$qform <- tmat
		hd$sform <- tmat

		quat1 <- matrixToQuatern(tmat)
		hd$quaternion <- quat1$quaternion
		hd$qfac <- quat1$qfac
		hd$pixdim[1] <- hd$qfac
		hd

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

	header$description <- character(80)
	header$auxfile <- character(24)

	header$qform_code <- 1
	header$sform_code <- 1
	header$quaternion <- NULL

	header$qoffset <- NULL
	header$qform <- NULL

	header$sform <- NULL
	header$intent_name <- character(16)
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


#' @keywords internal
#' @noRd
read_nifti_header <- function(fname) {
	header <- list()
	header$file_type <- "NIfTI"
	header$encoding <- "binary"
	header$version <- "1"

	conn <- NULL

	if (.isExtension(fname, ".nii") || .isExtension(fname, ".hdr")) {
		conn <- file(fname, open="rb")
	} else if (.isExtension(fname, ".nii.gz")) {
		conn <- gzfile(fname, open="rb")
		header$encoding <- "gzip"
	} else {
		stop(paste("illegal NIFTI header name", fname))
	}


	endian <- .getEndian(conn)

	header$file_name <- fname
	header$endian <- endian

	readBin(conn, what=integer(), n=10+18+4+2+1, size=1)

	#browser()

	header$diminfo <- readBin(conn, what=integer(), n=1, size=1)
	header$dimensions <- readBin(conn, integer(), n=8, size=2, endian=endian)

	header$dimensions[header$dimensions == 0] <- 1

	.checkDimensions(header$dimensions)

	header$num_dimensions <- header$dimensions[1]

	if (header$num_dimensions == 5 && header$dimensions[5] == 1) {
	  warning("read_nifti_header: collapsing fourth dimension because it has length = 1")
	  newdim <- header$dimensions[c(1:4, 6)]
	  header$num_dimensions <- 4
	  header$dimensions <- c(newdim, rep(1,3))
	  header$dimensions[1] <- 4
	}

	header$intent1 <-  readBin(conn, double(), n=1, size=4, endian=endian)
	header$intent2 <-  readBin(conn, double(), n=1, size=4, endian=endian)
	header$intent3 <-  readBin(conn, double(), n=1, size=4, endian=endian)


	header$intent_code <-  readBin(conn, integer(), n=1, size=2, endian=endian)
	header$datatype <- readBin(conn, integer(), n=1, size=2, endian=endian)
	header$data_storage <- .getDataStorage(header$datatype)
	header$bitpix <- readBin(conn, integer(), n=1, size=2, endian=endian)
	header$slice_start <- readBin(conn, integer(), n=1, size=2, endian=endian)
	header$pixdim <-  readBin(conn, double(), n=8, size=4, endian=endian)

	header$qfac = header$pixdim[1]

	if (header$qfac == 0) {
		header$qfac = 1
	}

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

	header$description <- readBin(conn, integer(), n=80, size=1, endian=endian)
	header$auxfile <- readBin(conn, integer(), n=24, size=1, endian=endian)

	header$qform_code <- readBin(conn, integer(), n=1, size=2, endian=endian)
	header$sform_code <- readBin(conn, integer(), n=1, size=2, endian=endian)
	header$quaternion <- readBin(conn, double(), n=3, size=4, endian=endian)

	header$qoffset <- readBin(conn, double(), n=3, size=4, endian=endian)
	header$qform <- quaternToMatrix(header$quaternion, header$qoffset, header$pixdim[2:4], header$qfac)

	sform  <- readBin(conn, double(), n=12, size=4, endian=endian)
	header$sform <- rbind(matrix(sform,3,4, byrow=T), c(0,0,0,1))
	header$intent_name <- readBin(conn, character(), n=16, size=1, endian=endian)
	header$magic <- suppressWarnings(readChar(conn, nchars=4))

	header$onefile <- F
	if (substr(header$magic,2,2) == "+") {
		header$onefile <- T
	}

	header$version <- substr(header$magic,3,3)

	close(conn)

	header

}


#' @keywords internal
#' @noRd
write_nifti_header <- function(niftiInfo, conn, close=TRUE) {
	endian <- niftiInfo$endian


	writeBin(as.integer(348), conn, 4, endian)
	writeBin(integer(34),conn,1,endian)
	writeChar("r", conn,1,eos=NULL)
	writeBin(as.integer(niftiInfo$diminfo), conn, size=1, endian) #diminfo, not supported currently -- write zero
	#writeBin(as.integer(niftiInfo$num_dimensions), conn, 2, endian)         #num dimensions

	stopifnot(length(niftiInfo$dimensions) == 8)
	stopifnot(length(niftiInfo$pixdim) == 8)

	writeBin(as.integer(niftiInfo$dimensions), conn, 2, endian)   #dimension vector
	writeBin(as.double(niftiInfo$intent1), conn, 4, endian)       #intent1
	writeBin(as.double(niftiInfo$intent2), conn, 4, endian)       #intent2
	writeBin(as.double(niftiInfo$intent3), conn, 4, endian)       #intent3
	writeBin(as.integer(niftiInfo$intent_code), conn, 2, endian)   #intent code
	writeBin(as.integer(niftiInfo$datatype),conn, 2, endian)      #datatype
	writeBin(as.integer(niftiInfo$bitpix),conn, 2, endian)        #bits per pixel
	writeBin(as.integer(niftiInfo$slice_start),conn, 2, endian)    #slice start
	writeBin(as.double(niftiInfo$pixdim), conn, 4, endian)        #pix dim
	writeBin(as.double(niftiInfo$vox_offset), conn, 4, endian)     #voxel offset
	writeBin(as.double(niftiInfo$scl_slope), conn, 4, endian)      #slope
	writeBin(as.double(niftiInfo$scl_intercept), conn, 4, endian)  #intercept
	writeBin(as.integer(niftiInfo$slice_end), conn, 2, endian)     #slice end
	writeBin(as.integer(niftiInfo$slice_code), conn, 1, endian)    #slice code
	writeBin(as.integer(niftiInfo$xyzt_units), conn, 1, endian)    #xyzt units
	writeBin(as.double(niftiInfo$cal_max), conn, 4, endian)        #cal max
	writeBin(as.double(niftiInfo$cal_min), conn, 4, endian)        #cal min
	writeBin(as.double(niftiInfo$slice_duration), conn, 4, endian) #slice duration
	writeBin(as.double(niftiInfo$toffset), conn, 4, endian)       #t offset
	writeBin(as.integer(niftiInfo$glmax), conn, 4, endian)        #glmax
	writeBin(as.integer(niftiInfo$glmin), conn, 4, endian)        #glmin
	writeBin(as.integer(niftiInfo$description), conn, 1, endian)  #description
	writeBin(as.integer(niftiInfo$auxfile), conn, 1, endian)      #aux_file
	writeBin(as.integer(niftiInfo$qform_code), conn, 2, endian)    #qform code
	writeBin(as.integer(niftiInfo$sform_code), conn, 2, endian)    #sform code

	writeBin(as.double(niftiInfo$quaternion), conn, 4, endian)    #quaternion
	writeBin(as.double(niftiInfo$qoffset), conn, 4, endian)       #qoffset
	writeBin(as.double(t(niftiInfo$sform[1:3,])), conn, 4, endian) #sform
	writeBin(as.integer(niftiInfo$intent_name), conn, 1, endian)    #intent_name
	writeChar(niftiInfo$magic, conn)                               #magic

	loc <- seek(conn)
	offset <- niftiInfo$vox_offset

	nbytes <- offset-loc

	## doesn't support extensions yet
	if (nbytes > 0) {
		writeBin(integer(nbytes), conn, size=1, endian)
	}

	if (close) {
		close(conn)
	}

	conn
}
