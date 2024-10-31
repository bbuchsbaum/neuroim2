#' @include all_class.R
NULL
#' @include all_generic.R
NULL


#' Create a Data Reader
#'
#' @description
#' This generic function creates a data reader for neuroimaging data. It is designed
#' to be flexible and work with various types of neuroimaging data formats.
#'
#' @param x An object specifying the information required to produce the reader.
#'   This could be a file path, a metadata object, or any other relevant information
#'   depending on the specific method implementation.
#' @param offset A numeric value specifying the byte offset (number of bytes to skip
#'   before reading). This allows for flexibility in reading data from different
#'   parts of a file or data stream.
#'
#' @return A data reader object. The exact type and structure of this object will
#'   depend on the specific method implementation.
#'
#' @details
#' The `data_reader` function is a generic that should be implemented for different
#' types of neuroimaging data formats (e.g., NIFTI, AFNI). Each implementation should
#' return an appropriate reader object that can be used to access the data.
#'
#' @seealso 
#' \code{\link{read_vol}}, \code{\link{FileMetaInfo-class}}
#'
#' @examples
#' \dontrun{
#' # Example usage (specific implementation may vary):
#' meta_info <- read_header("brain_scan.nii")
#' reader <- data_reader(meta_info, offset = 0)
#' }
#'
#' @export
#' @rdname data_reader-methods
setGeneric(name = "data_reader", 
           def = function(x, offset) standardGeneric("data_reader"))

#' dim of \code{FileMetaInfo}
#' @param x the object
#' @export
setMethod(f="dim", signature=signature("FileMetaInfo"),
		def=function(x) {
			x@dims
		})


# @rdname loadData-methods
# setMethod(f="loadData", signature=signature(""))


#' @rdname data_reader-methods
setMethod(f="data_reader", signature=signature("NIFTIMetaInfo"),
		def=function(x, offset=0) {
			if (x@descriptor@data_encoding == "gzip") {
				BinaryReader(gzfile(x@data_file, "rb"), x@data_offset+offset, .getRStorage(x@data_type), x@bytes_per_element, x@endian, .isSigned(x@data_type))
			} else {
				BinaryReader(x@data_file, x@data_offset+offset, .getRStorage(x@data_type), x@bytes_per_element, x@endian, .isSigned(x@data_type))
			}
		})

#' @rdname data_reader-methods
setMethod(f="data_reader", signature=signature("AFNIMetaInfo"),
		def=function(x, offset=0) {
			if (x@descriptor@data_encoding == "gzip") {
				BinaryReader(gzfile(x@data_file, "rb"), x@data_offset+offset, .getRStorage(x@data_type), x@bytes_per_element, x@endian,.isSigned(x@data_type))
			} else {
				BinaryReader(x@data_file, x@data_offset+offset, .getRStorage(x@data_type), x@bytes_per_element, x@endian,.isSigned(x@data_type))
			}
		})


#' @rdname trans-methods
setMethod(f="trans", signature=signature("MetaInfo"),
		def=function(x) {
			D <- min(length(x@dims), 3)
			trans <- diag(c(x@spacing,1))
			trans[1:D,D+1] <- x@origin
			trans
		})

#' @rdname trans-methods
setMethod(f="trans", signature=signature("NIFTIMetaInfo"),
		def=function(x) {
			x@header$qform
		})

#' @noRd 
#' @keywords internal
niftiDim <- function(nifti_header) {
	dimarray <- nifti_header$dimensions
	lastidx <- min(which(dimarray == 1)) - 1
	dimarray[2:lastidx]
}



#' Create a MetaInfo Object
#'
#' @description
#' This function creates a MetaInfo object, which contains metadata information for a neuroimaging volume.
#'
#' @param Dim A numeric vector specifying the image dimensions.
#' @param spacing A numeric vector specifying the voxel dimensions.
#' @param origin A numeric vector specifying the coordinate origin. Default is a zero vector of length equal to `spacing`.
#' @param data_type A character string specifying the type of the data (e.g., "FLOAT"). Default is "FLOAT".
#' @param label A character string specifying the name(s) of images. Default is an empty string.
#' @param spatial_axes An object specifying the image axes for spatial dimensions (x,y,z). 
#'   Default is OrientationList3D$AXIAL_LPI.
#' @param additional_axes An object specifying axes for dimensions > 3 (e.g., time, color band, direction). 
#'   Default is NullAxis.
#'
#' @return An instance of class \code{\linkS4class{MetaInfo}}.
#'
#' @details
#' The MetaInfo object encapsulates essential metadata about a neuroimaging volume, 
#' including its dimensions, voxel spacing, origin, data type, and axis orientations.
#' This information is crucial for correctly interpreting and manipulating the image data.
#'
#' @examples
#' # Create a MetaInfo object for a 3D image
#' meta <- MetaInfo(Dim = c(64, 64, 32), 
#'                  spacing = c(3, 3, 4), 
#'                  origin = c(0, 0, 0),
#'                  data_type = "FLOAT",
#'                  label = "T1_weighted")
#'
#' @seealso 
#' \code{\link{NIFTIMetaInfo}}, \code{\link{OrientationList3D}}
#'
#' @export
#' @rdname MetaInfo-class
MetaInfo <- function(Dim, spacing, origin = rep(0, length(spacing)), 
                     data_type = "FLOAT", label = "", 
                     spatial_axes = OrientationList3D$AXIAL_LPI, 
                     additional_axes = NullAxis) {
  new("MetaInfo",
      dims = Dim,
      spacing = spacing,
      origin = origin,
      data_type = data_type,
      label = label,
      spatial_axes = spatial_axes,
      additional_axes = additional_axes)
}



#' Create a NIFTIMetaInfo Object
#'
#' @description
#' This function creates a NIFTIMetaInfo object, which contains metadata information 
#' specific to NIFTI format neuroimaging files.
#'
#' @param descriptor An instance of class \code{\linkS4class{NIFTIFormat}}.
#' @param nifti_header A list returned by \code{readNIftiHeader} containing NIFTI header information.
#'
#' @return An instance of class \code{\linkS4class{NIFTIMetaInfo}}.
#'
#' @details
#' The NIFTIMetaInfo object extends the basic MetaInfo with additional fields specific to the NIFTI format.
#' It includes information about file locations, endianness, data offsets, and NIFTI-specific 
#' transformations and scaling factors.
#'
#' This function performs several steps:
#' 1. Checks that the input header is valid for a NIFTI file.
#' 2. Extracts relevant information from the NIFTI header.
#' 3. Computes derived information (e.g., bytes per element, dimensions).
#' 4. Creates and returns a new NIFTIMetaInfo object.
#'
#' @examples
#' \dontrun{
#' # Assuming you have a NIFTIFormat object and a NIFTI header
#' nifti_format <- NIFTIFormat()
#' nifti_header <- readNIftiHeader("brain.nii")
#' meta_info <- NIFTIMetaInfo(nifti_format, nifti_header)
#' }
#'
#' @seealso 
#' \code{\link{MetaInfo}}, \code{\link{NIFTIFormat-class}}, \code{\link{readNIftiHeader}}
#'
#' @export
#' @rdname FileMetaInfo-class
NIFTIMetaInfo <- function(descriptor, nifti_header) {
  stopifnot(!is.null(nifti_header$file_type) || (nifti_header$file_type == "NIFTI"))
  
  new("NIFTIMetaInfo",
      header_file = header_file(descriptor, nifti_header$file_name),
      data_file = data_file(descriptor, nifti_header$file_name),
      descriptor = descriptor,
      endian = nifti_header$endian,
      data_offset = nifti_header$vox_offset,
      data_type = nifti_header$data_storage,
      bytes_per_element = as.integer(.getDataSize(nifti_header$data_storage)),
      dims = niftiDim(nifti_header),
      spatial_axes = .nearestAnatomy(nifti_header$qform),
      additional_axes = NullAxis,
      spacing = nifti_header$pixdim[2:4],
      origin = nifti_header$qoffset,
      label = strip_extension(descriptor, basename(nifti_header$file_name)),
      intercept = nifti_header$scl_intercept,
      slope = nifti_header$scl_slope,
      header = nifti_header)
}





#' Create an AFNIMetaInfo Object
#'
#' @description
#' This function constructs an AFNIMetaInfo object, which contains metadata information 
#' specific to AFNI format neuroimaging files.
#'
#' @param descriptor An instance of class \code{\linkS4class{AFNIFormat}}.
#' @param afni_header A list returned by \code{read_afni_header} containing AFNI header information.
#'
#' @return An instance of class \code{\linkS4class{AFNIMetaInfo}}.
#'
#' @details
#' The AFNIMetaInfo object extends the basic MetaInfo with additional fields specific to the AFNI format.
#' It includes information about file locations, endianness, data types, dimensions, and AFNI-specific 
#' transformations and scaling factors.
#'
#' This function performs several steps:
#' 1. Extracts and processes dimension information from the AFNI header.
#' 2. Generates labels for each sub-brick if not provided in the header.
#' 3. Computes the transformation matrix from AFNI's IJK space to NIFTI's LPI space.
#' 4. Determines various metadata fields such as endianness, data type, and scaling factors.
#' 5. Creates and returns a new AFNIMetaInfo object.
#'
#' Note: The 'additional_axes' field is currently set to NullAxis, which may be incorrect 
#' for some AFNI datasets with more than 3 dimensions.
#'
#' @examples
#' \dontrun{
#' # Assuming you have an AFNIFormat object and an AFNI header
#' afni_format <- AFNIFormat()
#' afni_header <- read_afni_header("brain.HEAD")
#' meta_info <- AFNIMetaInfo(afni_format, afni_header)
#' }
#'
#' @seealso 
#' \code{\link{MetaInfo}}, \code{\link{AFNIFormat-class}}, \code{\link{read_afni_header}},
#' \code{\link{NIFTIMetaInfo}}
#'
#' @export
#' @rdname FileMetaInfo-class
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
		TLPI <- perm_mat(OrientationList3D$AXIAL_RAI) %*% Tdicom[1:3,]
    	TLPI <- rbind(TLPI, c(0,0,0,1))

		new("AFNIMetaInfo",
			header_file=header_file(descriptor, afni_header$file_name),
			data_file=data_file(descriptor, afni_header$file_name),
			descriptor=descriptor,
			endian=ifelse(afni_header[["BYTEORDER_STRING"]]$content == "MSB_FIRST", "big", "little"),
			data_offset=0,
			data_type=switch(afni_header$BRICK_TYPES$content[1], "0"="BYTE", "1"="SHORT", "3"="FLOAT"),
			bytes_per_element=as.integer(switch(afni_header$BRICK_TYPES$content[1], "0"=1, "1"=2, "3"=4)),
			dims=.Dim,
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

#' @export
setAs(from="MetaInfo", to="NIFTIMetaInfo", def=function(from) {
			if (inherits(from, "NIFTIMetaInfo")) {
				from
			} else {
				hdr <- as_nifti_header(from)
				desc <- find_descriptor(hdr$file_name)
				NIFTIMetaInfo(desc, hdr)
			}

		})


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
			cat("dimensions:", "\t", object@dims, "\n")
			cat("voxel size:", "\t", object@spacing, "\n")
			cat("origin:", "\t", object@origin, "\n")
			cat("label(s):", "\t", object@label, "\n")
			cat("intercept:", "\t", object@intercept, "\n")
			cat("slope:", "\t\t", object@slope, "\n\n")

			cat("additional format-specific info may be contained in @header slot", "\n")
		})




