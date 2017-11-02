#' @include all_class.R
NULL
#' @include all_generic.R
NULL


#' Generic function to create data reader
#' @param x an object specifying the infromation required to produce the reader
#' @param offset the byte offset (number of bytes to skip before reading)
#' @export data_reader
#' @rdname data_reader-methods
setGeneric(name="data_reader", def=function(x, offset) standardGeneric("data_reader"))

#' dim of \code{FileMetaInfo}
#' @param x the object
#' @export
setMethod(f="dim", signature=signature("FileMetaInfo"),
		def=function(x) {
			x@Dim
		})


# @rdname loadData-methods
# setMethod(f="loadData", signature=signature(""))


#' @rdname data_reader-methods
setMethod(f="data_reader", signature=signature("NIFTIMetaInfo"),
		def=function(x, offset=0) {
			if (x@fileDescriptor@dataEncoding == "gzip") {
				BinaryReader(gzfile(x@data_file, "rb"), x@data_offset+offset, .getRStorage(x@data_type), x@bytes_per_element, x@endian)
			} else {
				BinaryReader(x@data_file, x@data_offset+offset, .getRStorage(x@data_type), x@bytes_per_element, x@endian)
			}
		})

#' @rdname data_reader-methods
setMethod(f="data_reader", signature=signature("AFNIMetaInfo"),
		def=function(x, offset=0) {
			if (x@fileDescriptor@dataEncoding == "gzip") {
				BinaryReader(gzfile(x@data_file, "rb"), x@data_offset+offset, .getRStorage(x@data_type), x@bytes_per_element, x@endian)
			} else {
				BinaryReader(x@data_file, x@data_offset+offset, .getRStorage(x@data_type), x@bytes_per_element, x@endian)
			}
		})


#' @rdname trans-methods
setMethod(f="trans", signature=signature("MetaInfo"),
		def=function(x) {
			D <- min(length(x@Dim), 3)
			trans <- diag(c(x@spacing,1))
			trans[1:D,D+1] <- x@origin
			trans
		})

#' @rdname trans-methods
setMethod(f="trans", signature=signature("NIFTIMetaInfo"),
		def=function(x) {
			x@header$qform
		})

niftiDim <- function(nifti_header) {
	dimarray <- nifti_header$dimensions
	lastidx <- min(which(dimarray == 1)) - 1
	dimarray[2:lastidx]
}



#' This class contains meta information for an image
#'
#' @param Dim image dimensions
#' @param spacing voxel dimensions
#' @param origin coordinate origin
#' @param data_type the type of the data (e.g. "FLOAT")
#' @param label name(s) of images
#' @param spatial_axes image axes for spatial dimensions (x,y,z)
#' @param additional_axes axes for dimensions > 3 (e.g. time, color band, direction)
#' @return an instance of class \code{\linkS4class{MetaInfo}}
#' @export MetaInfo
#' @rdname MetaInfo-class
MetaInfo <- function(Dim, spacing, origin=rep(0, length(spacing)), data_type="FLOAT", label="", spatial_axes=OrientationList3D$AXIAL_LPI, additional_axes=NullAxis) {
	new("MetaInfo",
			Dim=Dim,
			spacing=spacing,
			origin=origin,
			data_type=data_type,
			label=label,
			spatial_axes=spatial_axes,
			additional_axes=additional_axes)
}



#' Constructor for \code{\linkS4class{NIFTIMetaInfo}} class
#' @param descriptor an instance of class \code{\linkS4class{NIFTIFormat}}
#' @param nifti_header a \code{list} returned by \code{readNIftiHeader}
#' @return an instance of class \code{\linkS4class{NIFTIMetaInfo}}
#' @export NIFTIMetaInfo
#' @rdname NIFTIMetaInfo-class
NIFTIMetaInfo <- function(descriptor, nifti_header) {
	stopifnot(!is.null(nifti_header$file_type) || (nifti_header$file_type == "NIFTI"))


	new("NIFTIMetaInfo",
			header_file=header_file(descriptor, nifti_header$file_name),
			data_file=data_file(descriptor, nifti_header$file_name),
			fileDescriptor=descriptor,
			endian=nifti_header$endian,
			data_offset=nifti_header$vox_offset,
			data_type=nifti_header$data_storage,
			bytes_per_element=as.integer(.getDataSize(nifti_header$data_storage)),
			Dim=niftiDim(nifti_header),
			spatial_axes=.nearestAnatomy(nifti_header$qform),
			additional_axes=NullAxis,
			spacing=nifti_header$pixdim[2:4],
			origin=nifti_header$qoffset,
			label=strip_extension(descriptor, basename(nifti_header$file_name)),
			intercept=nifti_header$sclIntercept,
			slope=nifti_header$sclSlope,
			header=nifti_header)
}



#' show a \code{FileMetaInfo}
#' @param object the object
#' @export
setMethod(f="show", signature=signature("FileMetaInfo"),
		def=function(object) {
			cat("an instance of class",  class(object), "\n\n")
			cat("header_file:", "\t", object@header_file, "\n")
			cat("data_file:", "\t", object@data_file, "\n")
			cat("endian:", "\t", object@endian, "\n")
			cat("data_offset:", "\t", object@data_offset, "\n")
			cat("data_type:", "\t", object@data_type, "\n")
			cat("dimensions:", "\t", object@Dim, "\n")
			cat("voxel size:", "\t", object@spacing, "\n")
			cat("origin:", "\t", object@origin, "\n")
			cat("label(s):", "\t", object@label, "\n")
			cat("intercept:", "\t", object@intercept, "\n")
			cat("slope:", "\t\t", object@slope, "\n\n")

			cat("additional format-specific info may be contained in @header slot", "\n")
		})


#' AFNIMetaInfo
#'
#' Constructor for \code{\linkS4class{AFNIMetaInfo}} class
#' @param descriptor an instance of class \code{\linkS4class{AFNIFormat}}
#' @param afni_header a \code{list} returned by \code{read_afni_header}
#' @return an instance of class \code{\linkS4class{AFNIMetaInfo}}
#' @export AFNIMetaInfo
#' @rdname AFNIMetaInfo-class
AFNIMetaInfo <- function(descriptor, afni_header) {
		.Dim <- afni_header$DATASET_DIMENSIONS$content[afni_header$DATASET_DIMENSIONS$content > 0]
		if (afni_header$DATASET_RANK$content[2] > 1) {
			.Dim <- c(.Dim, afni_header$DATASET_RANK$content[2])
		}


		labs <- if (is.null(afni_header$BRICK_LABS$content)) {
			labs <- paste("#", seq(0, afni_header$DATASET_RANK$content[2]), sep="")
		} else {
			afni_header$BRICK_LABS$content
		}


    ## AFNI contains a transform from IJK to dicom (RAI) space.
    ## We want the transform to go from IJK to nifti (LPI) space
		Tdicom <- matrix(afni_header$IJK_TO_DICOM$content, 3,4, byrow=TRUE)
		TLPI <- permMat(OrientationList3D$AXIAL_RAI) %*% Tdicom[1:3,]
    TLPI <- rbind(TLPI, c(0,0,0,1))

		new("AFNIMetaInfo",
			header_file=header_file(descriptor, afni_header$file_name),
			data_file=data_file(descriptor, afni_header$file_name),
			fileDescriptor=descriptor,
			endian=ifelse(afni_header[["BYTEORDER_STRING"]]$content == "MSB_FIRST", "big", "little"),
			data_offset=0,
			data_type=switch(afni_header$BRICK_TYPES$content[1], "0"="BYTE", "1"="SHORT", "3"="FLOAT"),
			bytes_per_element=as.integer(switch(afni_header$BRICK_TYPES$content[1], "0"=1, "1"=2, "3"=4)),
			Dim=.Dim,
			spatial_axes=.nearestAnatomy(TLPI),
			additional_axes=NullAxis,            # incorrect
			spacing=abs(afni_header$DELTA$content),
			origin=afni_header$ORIGIN$content,
			label=labs,
			intercept=0,
			slope=ifelse(afni_header$BRICK_FLOAT_FACS$content == 0, 1, afni_header$BRICK_FLOAT_FACS$content),
			header=afni_header)
}


#' read header information of an image file
#'
#'
#' @param file_name the name of the file to read
#' @return an instance of class \code{\linkS4class{FileMetaInfo}}
#' @export read_header
read_header <- function(file_name) {
	desc <- find_descriptor(file_name)
	if (is.null(desc)) {
		stop(paste("could not find reader for file: ", file_name))
	}

	read_meta_info(desc, file_name)
}

setAs(from="MetaInfo", to="NIFTIMetaInfo", def=function(from) {
			if (inherits(from, "NIFTIMetaInfo")) {
				from
			} else {
				hdr <- as.nifti.header(from)
				desc <- find_descriptor(hdr$file_name)
				NIFTIMetaInfo(desc, hdr)
			}

		})



