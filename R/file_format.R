#' @include all_class.R
NULL
#' @include all_generic.R


#' Check if a File Matches the FileFormat
#'
#' @description
#' This method checks if a given file name matches the specified FileFormat,
#' verifying the existence of both header and data files.
#'
#' @param x A FileFormat object.
#' @param file_name A character string representing the file name to be matched.
#'
#' @return A logical value: TRUE if the file matches the format and both header 
#'   and data files exist, FALSE otherwise.
#'
#' @details
#' The function performs the following checks:
#' 1. If the file name matches the header file format, it checks for the existence of the corresponding data file.
#' 2. If the file name matches the data file format, it checks for the existence of the corresponding header file.
#' 3. If neither condition is met, it returns FALSE.
#'
#' @seealso 
#' \code{\link{FileFormat-class}}, \code{\link{header_file_matches}}, \code{\link{data_file_matches}}
#'
#' @examples
#' \dontrun{
#' file_format <- new("FileFormat", header_extension = "hdr", data_extension = "img")
#' file_matches(file_format, "brain_scan.hdr")
#' file_matches(file_format, "brain_scan.img")
#' }
#'
#' @export
#' @rdname file_matches-methods
setMethod(f = "file_matches", 
          signature = signature(x = "FileFormat", file_name = "character"),
          def = function(x, file_name) {
            if (header_file_matches(x, file_name)) {
              file.exists(paste(strip_extension(x, file_name), ".", x@data_extension, sep = ""))
            } else if (data_file_matches(x, file_name)) {
              file.exists(paste(strip_extension(x, file_name), ".", x@header_extension, sep = ""))
            } else {
              FALSE
            }
          })

#' Check if a File Matches the Header File Format
#'
#' @description
#' This method checks if a given file name conforms to the header file format
#' specified in the FileFormat object.
#'
#' @param x A FileFormat object.
#' @param file_name A character string representing the file name to be matched.
#'
#' @return A logical value: TRUE if the file name matches the header file format, 
#'   FALSE otherwise.
#'
#' @details
#' The function uses a regular expression to check if the file name ends with 
#' the header extension specified in the FileFormat object. It returns TRUE if 
#' there's a match, and FALSE otherwise.
#'
#' @seealso 
#' \code{\link{FileFormat-class}}, \code{\link{file_matches}}, \code{\link{data_file_matches}}
#'
#' @examples
#' \dontrun{
#' file_format <- new("FileFormat", header_extension = "hdr", data_extension = "img")
#' header_file_matches(file_format, "brain_scan.hdr")
#' header_file_matches(file_format, "brain_scan.img")
#' }
#'
#' @export
#' @rdname header_file_matches-methods
setMethod(f = "header_file_matches", 
          signature = signature(x = "FileFormat", file_name = "character"),
          def = function(x, file_name) {
            regexpr(paste(".*", x@header_extension, "$", sep = ""), file_name) > 0
          })

#' Check if a File Matches the Data File Format
#'
#' @description
#' This method checks if a given file name conforms to the data file format
#' specified in the FileFormat object.
#'
#' @param x A FileFormat object.
#' @param file_name A character string representing the file name to be matched.
#'
#' @return A logical value: TRUE if the file name matches the data file format, 
#'   FALSE otherwise.
#'
#' @details
#' The function uses a regular expression to check if the file name ends with 
#' the data extension specified in the FileFormat object. It returns TRUE if 
#' there's a match, and FALSE otherwise.
#'
#' @seealso 
#' \code{\link{FileFormat-class}}, \code{\link{file_matches}}, \code{\link{header_file_matches}}
#'
#' @examples
#' \dontrun{
#' file_format <- new("FileFormat", header_extension = "hdr", data_extension = "img")
#' data_file_matches(file_format, "brain_scan.img")
#' data_file_matches(file_format, "brain_scan.hdr")
#' }
#'
#' @export
#' @rdname data_file_matches-methods
setMethod(f = "data_file_matches", 
          signature = signature(x = "FileFormat", file_name = "character"),
          def = function(x, file_name) {
            regexpr(paste(".*", x@data_extension, "$", sep = ""), file_name) > 0
          })

#' Get the Header File Name
#'
#' @description
#' This method derives the header file name from a given file name based on 
#' the FileFormat specifications.
#'
#' @param x A FileFormat object.
#' @param file_name A character string representing the file name to be processed.
#'
#' @return A character string representing the correct header file name.
#'
#' @details
#' The function performs the following steps:
#' 1. If the input file_name already matches the header file format, it returns the file_name as is.
#' 2. If the file_name matches the data file format, it constructs and returns the corresponding header file name.
#' 3. If the file_name doesn't match either format, it throws an error.
#'
#' @seealso 
#' \code{\link{FileFormat-class}}, \code{\link{header_file_matches}}, \code{\link{data_file_matches}}, \code{\link{strip_extension}}
#'
#' @examples
#' \dontrun{
#' file_format <- new("FileFormat", header_extension = "hdr", data_extension = "img")
#' header_file(file_format, "brain_scan.hdr")
#' header_file(file_format, "brain_scan.img")
#' }
#'
#' @export
#' @rdname header_file-methods
setMethod(f = "header_file",
          signature = signature(x = "FileFormat", file_name = "character"),
          def = function(x, file_name) {
            if (header_file_matches(x, file_name)) {
              file_name
            } else if (data_file_matches(x, file_name)) {
              paste(strip_extension(x, file_name), x@header_extension, sep = ".")
            } else {
              stop(paste("Could not derive header file name from:", file_name))
            }
          })

#' Get the Data File Name
#'
#' @description
#' This method derives the data file name from a given file name based on 
#' the FileFormat specifications.
#'
#' @param x A FileFormat object.
#' @param file_name A character string representing the file name to be processed.
#'
#' @return A character string representing the correct data file name.
#'
#' @details
#' The function performs the following steps:
#' 1. If the input file_name already matches the data file format, it returns the file_name as is.
#' 2. If the file_name matches the header file format, it constructs and returns the corresponding data file name.
#' 3. If the file_name doesn't match either format, it throws an error.
#'
#' @seealso 
#' \code{\link{FileFormat-class}}, \code{\link{data_file_matches}}, \code{\link{header_file_matches}}, \code{\link{strip_extension}}
#'
#' @examples
#' \dontrun{
#' file_format <- new("FileFormat", header_extension = "hdr", data_extension = "img")
#' data_file(file_format, "brain_scan.img")
#' data_file(file_format, "brain_scan.hdr")
#' }
#'
#' @export
#' @rdname data_file-methods
setMethod(f = "data_file",
          signature = signature(x = "FileFormat", file_name = "character"),
          def = function(x, file_name) {
            if (data_file_matches(x, file_name)) {
              file_name
            } else if (header_file_matches(x, file_name)) {
              paste(strip_extension(x, file_name), x@data_extension, sep = ".")
            } else {
              stop(paste("Could not derive data file name from:", file_name))
            }
          })

#' Strip File Extension
#'
#' @description
#' This method removes the file extension from a given file name based on 
#' the FileFormat specifications.
#'
#' @param x A FileFormat object.
#' @param file_name A character string representing the file name to be processed.
#'
#' @return A character string representing the file name without the extension.
#'
#' @details
#' The function performs the following steps:
#' 1. If the file_name matches the header file format, it removes the header extension.
#' 2. If the file_name matches the data file format, it removes the data extension.
#' 3. If the file_name doesn't match either format, it throws an error.
#'
#' @seealso 
#' \code{\link{FileFormat-class}}, \code{\link{header_file_matches}}, \code{\link{data_file_matches}}
#'
#' @examples
#' \dontrun{
#' file_format <- new("FileFormat", header_extension = "hdr", data_extension = "img")
#' strip_extension(file_format, "brain_scan.hdr")
#' strip_extension(file_format, "brain_scan.img")
#' }
#'
#' @export
#' @rdname strip_extension-methods
setMethod(f = "strip_extension",
          signature = signature(x = "FileFormat", file_name = "character"),
          def = function(x, file_name) {
            if (header_file_matches(x, file_name)) {
              ret <- strsplit(file_name, paste(x@header_extension, "$", sep = ""))[[1]][1]
              substr(ret, 1, nchar(ret) - 1)
            } else if (data_file_matches(x, file_name)) {
              ret <- strsplit(file_name, paste(x@data_extension, "$", sep = ""))[[1]][1]
              substr(ret, 1, nchar(ret) - 1)
            } else {
              stop(paste("File does not match descriptor:", x))
            }
          })



#' @keywords internal
#' @noRd
.read_meta_info <- function(desc, file_name, read_func, constructor) {

	hfile <- header_file(desc, file_name)
	header <- read_func(hfile)
	header$file_name <- hfile
	constructor(desc, header)
}

#' Read Image Meta Information
#'
#' @description
#' These methods read meta information from image files based on their format (NIFTI or AFNI).
#'
#' @param x A FileFormat object (either NIFTIFormat or AFNIFormat).
#' @param file_name A character string representing the file name to read meta information from.
#'
#' @return An object of class NIFTIMetaInfo or AFNIMetaInfo, depending on the input format.
#'
#' @details
#' These methods use format-specific functions to read the header information and 
#' create the appropriate meta information object. The `.read_meta_info` helper function 
#' is used internally to streamline the process for both formats.
#'
#' @seealso 
#' \code{\link{FileFormat-class}}, \code{\link{NIFTIFormat-class}}, \code{\link{AFNIFormat-class}},
#' \code{\link{NIFTIMetaInfo-class}}, \code{\link{AFNIMetaInfo-class}}
#'
#' @examples
#' \dontrun{
#' # For NIFTI format
#' nifti_format <- new("NIFTIFormat")
#' nifti_meta <- read_meta_info(nifti_format, "brain_scan.nii")
#'
#' # For AFNI format
#' afni_format <- new("AFNIFormat")
#' afni_meta <- read_meta_info(afni_format, "brain_scan+orig.HEAD")
#' }
#'
#' @export
#' @rdname read_meta_info-methods
setMethod(f="read_meta_info",signature=signature(x= "NIFTIFormat"),
		def=function(x, file_name) {
			.read_meta_info(x, file_name, read_nifti_header, NIFTIMetaInfo)
		})

#' @export
#' @rdname read_meta_info-methods
setMethod(f="read_meta_info",signature=signature(x= "AFNIFormat"),
		def=function(x, file_name) {
			.read_meta_info(x, file_name, read_afni_header, AFNIMetaInfo)

		})


#' @keywords internal
#' @noRd
find_descriptor <- function(file_name) {
	if (file_matches(NIFTI, file_name)) NIFTI
	else if (file_matches(NIFTI_GZ, file_name)) NIFTI_GZ
	else if (file_matches(NIFTI_PAIR, file_name)) NIFTI_PAIR
	else if (file_matches(NIFTI_PAIR_GZ, file_name)) NIFTI_PAIR_GZ
	else if (file_matches(AFNI, file_name)) AFNI
	else if (file_matches(AFNI_GZ, file_name)) AFNI_GZ
	else NULL
}


#' @keywords internal
#' @noRd
AFNI <- new("AFNIFormat",
		file_format="AFNI",
		header_encoding="raw",
		header_extension="HEAD",
		data_encoding="raw",
		data_extension="BRIK")


#' @keywords internal
#' @noRd
AFNI_GZ <- new("AFNIFormat",
		file_format="AFNI",
		header_encoding="gzip",
		header_extension="HEAD",
		data_encoding="gzip",
		data_extension="BRIK.gz")


#' @keywords internal
#' @noRd
NIFTI <- new("NIFTIFormat",
		file_format="NIFTI",
		header_encoding="raw",
		header_extension="nii",
		data_encoding="raw",
		data_extension="nii")


#' @keywords internal
#' @noRd
NIFTI_GZ <- new("NIFTIFormat",
		file_format="NIFTI",
		header_encoding="gzip",
		header_extension="nii.gz",
		data_encoding="gzip",
		data_extension="nii.gz")


#' @keywords internal
#' @noRd
NIFTI_PAIR <- new("NIFTIFormat",
		file_format="NIFTI",
		header_encoding="raw",
		header_extension="hdr",
		data_encoding="raw",
		data_extension="img")


#' @keywords internal
NIFTI_PAIR_GZ <- new("NIFTIFormat",
		file_format="NIFTI",
		header_encoding="gzip",
		header_extension="hdr.gz",
		data_encoding="gzip",
		data_extension="img.gz")


