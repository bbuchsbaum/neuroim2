#' @include all_class.R
#' @include all_generic.R
NULL

#' File Format Operations for Neuroimaging Data
#'
#' @name FileFormat-operations
#' @description
#' A collection of methods for handling neuroimaging file formats with separate header
#' and data files (e.g., ANALYZE, NIFTI). These methods provide functionality for
#' file name validation, extension handling, and file path manipulation.
#'
#'
#' @section File Format Structure:
#' Neuroimaging formats often use paired files:
#' \itemize{
#'   \item A header file (e.g., '.hdr') containing metadata
#'   \item A data file (e.g., '.img') containing the actual image data
#' }
#'
#' @section Common Operations:
#' \itemize{
#'   \item Validating file names against format specifications
#'   \item Converting between header and data file names
#'   \item Checking file existence and compatibility
#' }
#'
#' @importFrom methods new
#' @importFrom assertthat assert_that
#'
NULL

#' Check if a File Matches the FileFormat
#'
#' @description
#' Validates whether a file name conforms to the specified FileFormat and verifies
#' the existence of both header and data files.
#'
#' @param x A \linkS4class{FileFormat} object specifying the format requirements
#' @param file_name A character string specifying the file name to validate
#'
#' @return A logical value: \code{TRUE} if the file matches the format and both header
#'   and data files exist, \code{FALSE} otherwise
#'
#' @details
#' The function performs the following validation steps:
#' \enumerate{
#'   \item Checks if the file name matches either the header or data format
#'   \item Verifies the existence of the corresponding paired file
#'   \item Returns \code{FALSE} if either check fails
#' }
#'
#' File names are validated using case-sensitive extension matching.
#'
#' @examples
#' \dontrun{
#' # Create a FileFormat for ANALYZE format
#' fmt <- new("FileFormat", header_extension = "hdr", data_extension = "img")
#'
#' # Check if files exist and match format
#' file_matches(fmt, "brain_scan.hdr")  # TRUE if both .hdr and .img exist
#' file_matches(fmt, "brain_scan.nii")  # FALSE for wrong extension
#' }
#'
#' @seealso
#' \code{\link{header_file_matches}}, \code{\link{data_file_matches}} for individual
#' file type checking
#'
#' @export
#' @rdname file_matches-methods
setMethod(f = "file_matches",
          signature = signature(x = "FileFormat", file_name = "character"),
          def = function(x, file_name) {
            assert_that(is.character(file_name) && length(file_name) == 1,
                       msg = "'file_name' must be a single character string")
            assert_that(!is.na(file_name) && nchar(file_name) > 0,
                       msg = "'file_name' cannot be NA or empty")

            # Check file existence and format matching
            if (header_file_matches(x, file_name)) {
              data_file <- paste(strip_extension(x, file_name), ".", x@data_extension, sep = "")
              file.exists(data_file)
            } else if (data_file_matches(x, file_name)) {
              header_file <- paste(strip_extension(x, file_name), ".", x@header_extension, sep = "")
              file.exists(header_file)
            } else {
              FALSE
            }
          })



#' Check if a File Matches the Header File Format
#'
#' @description
#' Validates whether a file name conforms to the header file format specification.
#'
#' @param x A \linkS4class{FileFormat} object specifying the format requirements
#' @param file_name A character string specifying the file name to validate
#'
#' @return A logical value: \code{TRUE} if the file name matches the header format,
#'   \code{FALSE} otherwise
#'
#' @details
#' The function performs case-sensitive pattern matching to verify that the file name
#' ends with the specified header extension. The match is performed using a regular
#' expression that ensures the extension appears at the end of the file name.
#'
#' @examples
#' \dontrun{
#' fmt <- new("FileFormat", header_extension = "hdr", data_extension = "img")
#' header_file_matches(fmt, "brain_scan.hdr")  # TRUE
#' header_file_matches(fmt, "brain_scan.img")  # FALSE
#' header_file_matches(fmt, "brain.hdr.gz")    # FALSE
#' }
#'
#' @seealso
#' \code{\link{file_matches}}, \code{\link{data_file_matches}} for related
#' file format validation
#'
#' @export
#' @rdname header_file_matches-methods
setMethod(f = "header_file_matches",
          signature = signature(x = "FileFormat", file_name = "character"),
          def = function(x, file_name) {
            # Input validation
            if (!is.character(file_name) || length(file_name) != 1) {
              stop("'file_name' must be a single character string")
            }
            if (is.na(file_name) || nchar(file_name) == 0) {
              stop("'file_name' cannot be NA or empty")
            }

            # Use anchored regex for exact extension matching
            pattern <- paste0("\\.", x@header_extension, "$")
            grepl(pattern, file_name, perl = TRUE)
          })

#' Check if a File Matches the Data File Format
#'
#' @description
#' Validates whether a file name conforms to the data file format specification.
#'
#' @param x A \linkS4class{FileFormat} object specifying the format requirements
#' @param file_name A character string specifying the file name to validate
#'
#' @return A logical value: \code{TRUE} if the file name matches the data format,
#'   \code{FALSE} otherwise
#'
#' @details
#' The function performs case-sensitive pattern matching to verify that the file name
#' ends with the specified data extension. The match is performed using a regular
#' expression that ensures the extension appears at the end of the file name.
#'
#' @examples
#' \dontrun{
#' fmt <- new("FileFormat", header_extension = "hdr", data_extension = "img")
#' data_file_matches(fmt, "brain_scan.img")  # TRUE
#' data_file_matches(fmt, "brain_scan.hdr")  # FALSE
#' data_file_matches(fmt, "brain.img.gz")    # FALSE
#' }
#'
#' @seealso
#' \code{\link{file_matches}}, \code{\link{header_file_matches}} for related
#' file format validation
#'
#' @export
#' @rdname data_file_matches-methods
setMethod(f = "data_file_matches",
          signature = signature(x = "FileFormat", file_name = "character"),
          def = function(x, file_name) {
            # Input validation
            if (!is.character(file_name) || length(file_name) != 1) {
              stop("'file_name' must be a single character string")
            }
            if (is.na(file_name) || nchar(file_name) == 0) {
              stop("'file_name' cannot be NA or empty")
            }

            # Use anchored regex for exact extension matching
            pattern <- paste0("\\.", x@data_extension, "$")
            grepl(pattern, file_name, perl = TRUE)
          })

#' Get the Header File Name
#'
#' @description
#' Derives the header file name from a given file name based on the FileFormat
#' specifications.
#'
#' @param x A \linkS4class{FileFormat} object specifying the format requirements
#' @param file_name A character string specifying the file name to derive the header
#'   file name from
#'
#' @return A character string representing the header file name
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item If the input file_name already matches the header file format, it returns
#'     the file_name as is.
#'   \item If the file_name matches the data file format, it constructs and returns
#'     the corresponding header file name.
#'   \item If the file_name doesn't match either format, it throws an error.
#' }
#'
#' @examples
#' \dontrun{
#' fmt <- new("FileFormat", header_extension = "hdr", data_extension = "img")
#' header_file(fmt, "brain_scan.hdr")  # Returns "brain_scan.hdr"
#' header_file(fmt, "brain_scan.img")  # Returns "brain_scan.hdr"
#' }
#'
#' @seealso
#' \code{\link{data_file}}, \code{\link{strip_extension}} for related file name
#' manipulation
#'
#' @export
#' @rdname header_file-methods
setMethod(f = "header_file",
          signature = signature(x = "FileFormat", file_name = "character"),
          def = function(x, file_name) {
            # Input validation
            if (!is.character(file_name) || length(file_name) != 1) {
              stop("'file_name' must be a single character string")
            }
            if (is.na(file_name) || nchar(file_name) == 0) {
              stop("'file_name' cannot be NA or empty")
            }

            # Derive header file name based on format
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
#' Derives the data file name from a given file name based on the FileFormat
#' specifications.
#'
#' @param x A \linkS4class{FileFormat} object specifying the format requirements
#' @param file_name A character string specifying the file name to derive the data
#'   file name from
#'
#' @return A character string representing the data file name
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item If the input file_name already matches the data file format, it returns
#'     the file_name as is.
#'   \item If the file_name matches the header file format, it constructs and returns
#'     the corresponding data file name.
#'   \item If the file_name doesn't match either format, it throws an error.
#' }
#'
#' @examples
#' \dontrun{
#' fmt <- new("FileFormat", header_extension = "hdr", data_extension = "img")
#' data_file(fmt, "brain_scan.img")  # Returns "brain_scan.img"
#' data_file(fmt, "brain_scan.hdr")  # Returns "brain_scan.img"
#' }
#'
#' @seealso
#' \code{\link{header_file}}, \code{\link{strip_extension}} for related file name
#' manipulation
#'
#' @export
#' @rdname data_file-methods
setMethod(f = "data_file",
          signature = signature(x = "FileFormat", file_name = "character"),
          def = function(x, file_name) {
            # Input validation
            if (!is.character(file_name) || length(file_name) != 1) {
              stop("'file_name' must be a single character string")
            }
            if (is.na(file_name) || nchar(file_name) == 0) {
              stop("'file_name' cannot be NA or empty")
            }

            # Derive data file name based on format
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
#' Removes the file extension from a given file name based on the FileFormat
#' specifications.
#'
#' @param x A \linkS4class{FileFormat} object specifying the format requirements
#' @param file_name A character string specifying the file name to strip the
#'   extension from
#'
#' @return A character string representing the file name without the extension
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item If the file_name matches the header file format, it removes the header
#'     extension.
#'   \item If the file_name matches the data file format, it removes the data
#'     extension.
#'   \item If the file_name doesn't match either format, it throws an error.
#' }
#'
#' @examples
#' \dontrun{
#' fmt <- new("FileFormat", header_extension = "hdr", data_extension = "img")
#' strip_extension(fmt, "brain_scan.hdr")  # Returns "brain_scan"
#' strip_extension(fmt, "brain_scan.img")  # Returns "brain_scan"
#' }
#'
#' @seealso
#' \code{\link{header_file}}, \code{\link{data_file}} for related file name
#' manipulation
#'
#' @export
#' @rdname strip_extension-methods
setMethod(f = "strip_extension",
          signature = signature(x = "FileFormat", file_name = "character"),
          def = function(x, file_name) {
             # Input validation
            if (!is.character(file_name) || length(file_name) != 1) {
              stop("'file_name' must be a single character string")
            }
            if (is.na(file_name) || nchar(file_name) == 0) {
              stop("'file_name' cannot be NA or empty")
            }

            # Remove extension based on format
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
#' Reads meta information from image files based on their format (NIFTI or AFNI).
#'
#' @param x A \linkS4class{FileFormat} object (either NIFTIFormat or AFNIFormat)
#' @param file_name A character string specifying the file name to read meta
#'   information from
#'
#' @return An object of class \linkS4class{NIFTIMetaInfo} or \linkS4class{AFNIMetaInfo},
#'   depending on the input format
#'
#' @details
#' These methods use format-specific functions to read the header information and
#' create the appropriate meta information object. The `.read_meta_info` helper
#' function is used internally to streamline the process for both formats.
#'
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
setMethod(f = "read_meta_info", signature = signature(x = "NIFTIFormat"),
          def = function(x, file_name) {
            .read_meta_info(x, file_name, read_nifti_header, NIFTIMetaInfo)
          })

#' @export
#' @rdname read_meta_info-methods
setMethod(f = "read_meta_info", signature = signature(x = "AFNIFormat"),
          def = function(x, file_name) {
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
            file_format = "AFNI",
            header_encoding = "raw",
            header_extension = "HEAD",
            data_encoding = "raw",
            data_extension = "BRIK")

#' @keywords internal
#' @noRd
AFNI_GZ <- new("AFNIFormat",
               file_format = "AFNI",
               header_encoding = "gzip",
               header_extension = "HEAD",
               data_encoding = "gzip",
               data_extension = "BRIK.gz")

#' @keywords internal
#' @noRd
NIFTI <- new("NIFTIFormat",
             file_format = "NIFTI",
             header_encoding = "raw",
             header_extension = "nii",
             data_encoding = "raw",
             data_extension = "nii")

#' @keywords internal
#' @noRd
NIFTI_GZ <- new("NIFTIFormat",
                file_format = "NIFTI",
                header_encoding = "gzip",
                header_extension = "nii.gz",
                data_encoding = "gzip",
                data_extension = "nii.gz")

#' @keywords internal
#' @noRd
NIFTI_PAIR <- new("NIFTIFormat",
                  file_format = "NIFTI",
                  header_encoding = "raw",
                  header_extension = "hdr",
                  data_encoding = "raw",
                  data_extension = "img")

#' @keywords internal
NIFTI_PAIR_GZ <- new("NIFTIFormat",
                     file_format = "NIFTI",
                     header_encoding = "gzip",
                     header_extension = "hdr.gz",
                     data_encoding = "gzip",
                     data_extension = "img.gz")
