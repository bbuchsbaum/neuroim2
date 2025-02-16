#' @include all_class.R
NULL
#' @include all_generic.R
NULL

#' Create a Data Reader for Neuroimaging Data
#'
#' @title Create a Data Reader
#' @description Creates a data reader for accessing neuroimaging data from various file formats.
#'   The reader provides a unified interface for reading data regardless of the underlying format.
#'
#' @param x An object containing metadata required to create the reader (e.g., file path, format info)
#' @param offset Numeric. Byte offset where data reading should begin. Default is 0.
#'
#' @return A BinaryReader object configured for the specific data format
#'
#' @details
#' The data_reader function is a generic that creates appropriate readers for different
#' neuroimaging formats. It handles:
#' \itemize{
#'   \item File format detection and validation
#'   \item Endianness configuration
#'   \item Data type conversion
#'   \item Compression handling (e.g., gzip)
#'   \item Proper byte alignment
#' }
#'
#' @examples
#' \dontrun{
#' # Create reader for NIFTI file
#' meta <- read_header("brain.nii")
#' reader <- data_reader(meta, offset = 0)
#'
#' # Read first 100 voxels
#' data <- reader$read(100)
#' }
#'
#' @seealso
#' \code{\link{read_header}} for reading headers,
#' \code{\linkS4class{BinaryReader}} for reading binary data
#'
#' @importFrom methods new setGeneric
#'
#' @export
setGeneric("data_reader", function(x, offset) standardGeneric("data_reader"))

#' Get Dimensions of FileMetaInfo Object
#'
#' @param x A FileMetaInfo object
#' @return Integer vector of image dimensions
#' @rdname dim-methods
#' @export
setMethod("dim", "FileMetaInfo", function(x) x@dims)


#' @param offset the offset to read from the file
#' @rdname data_reader-methods
setMethod("data_reader", "NIFTIMetaInfo",
  function(x, offset = 0) {
    assert_that(is.numeric(offset) && length(offset) == 1,
                msg = "'offset' must be a single numeric value")

    total_offset <- x@data_offset + offset

    if (x@descriptor@data_encoding == "gzip") {
      con <- tryCatch({
        gzfile(x@data_file, "rb")
      }, error = function(e) {
        stop("Failed to open gzipped file: ", x@data_file, "\n", e$message)
      })
    } else {
      con <- x@data_file
    }

    BinaryReader(con,
                total_offset,
                .getRStorage(x@data_type),
                x@bytes_per_element,
                x@endian,
                .isSigned(x@data_type))
  })

#' Create Data Reader for AFNI Format
#'
#' @param x AFNIMetaInfo object
#' @param offset Numeric byte offset
#' @return BinaryReader object
#' @rdname data_reader-methods
setMethod("data_reader", "AFNIMetaInfo",
  function(x, offset = 0) {
    assert_that(is.numeric(offset) && length(offset) == 1,
                msg = "'offset' must be a single numeric value")

    total_offset <- x@data_offset + offset

    if (x@descriptor@data_encoding == "gzip") {
      con <- tryCatch({
        gzfile(x@data_file, "rb")
      }, error = function(e) {
        stop("Failed to open gzipped file: ", x@data_file, "\n", e$message)
      })
    } else {
      con <- x@data_file
    }

    BinaryReader(con,
                total_offset,
                .getRStorage(x@data_type),
                x@bytes_per_element,
                x@endian,
                .isSigned(x@data_type))
  })

#' Get transformation matrix
#' @name trans
#' @rdname trans-methods
#' @aliases trans,MetaInfo-method trans,NIFTIMetaInfo-method
#' @export
setMethod("trans", "MetaInfo",
  function(x) {
    D <- min(length(x@dims), 3)
    trans <- diag(c(x@spacing, 1))
    trans[1:D, D + 1] <- x@origin
    trans
  })

#' Get NIFTI Transformation Matrix
#'
#' @param x NIFTIMetaInfo object
#' @return 4x4 transformation matrix
#' @keywords internal
#' @noRd
setMethod("trans", "NIFTIMetaInfo",
  function(x) x@header$qform)

#' Extract NIFTI Dimensions
#'
#' @param nifti_header NIFTI header list
#' @return Vector of image dimensions
#' @keywords internal
#' @noRd
niftiDim <- function(nifti_header) {
  dimarray <- nifti_header$dimensions
  assert_that(is.numeric(dimarray),
              msg = "Invalid dimension array in NIFTI header")
  lastidx <- min(which(dimarray == 1)) - 1
  assert_that(lastidx >= 1,
              msg = "Invalid dimension specification in NIFTI header")
  dimarray[2:lastidx]
}

#' Create MetaInfo Object
#'
#' @title Create Neuroimaging Metadata Object
#' @description Creates a MetaInfo object containing essential metadata for neuroimaging data,
#'   including dimensions, spacing, orientation, and data type information.
#'
#' @param Dim Integer vector. Image dimensions (e.g., c(64, 64, 32) for 3D).
#' @param spacing Numeric vector. Voxel dimensions in mm.
#' @param origin Numeric vector. Coordinate origin. Default is zero vector.
#' @param data_type Character. Data type (e.g., "FLOAT", "SHORT"). Default is "FLOAT".
#' @param label Character. Image label(s). Default is "".
#' @param spatial_axes Object. Spatial orientation. Default is OrientationList3D$AXIAL_LPI.
#' @param additional_axes Object. Non-spatial axes. Default is NullAxis.
#'
#' @return A MetaInfo object
#'
#' @details
#' The MetaInfo object is fundamental for:
#' \itemize{
#'   \item Spatial interpretation of image data
#'   \item Data type handling and conversion
#'   \item Memory allocation and mapping
#'   \item File I/O operations
#' }
#'
#' Input validation ensures:
#' \itemize{
#'   \item Dimensions are positive integers
#'   \item Spacing values are positive
#'   \item Origin coordinates are finite
#'   \item Data type is supported
#' }
#'
#' @examples
#' # Create metadata for 3D structural MRI
#' meta <- MetaInfo(
#'   Dim = c(256, 256, 180),
#'   spacing = c(1, 1, 1),
#'   data_type = "FLOAT",
#'   label = "T1w"
#' )
#'
#' # Get image dimensions
#' dim(meta)
#'
#' # Get transformation matrix
#' trans(meta)
#'
#' @seealso
#' \code{\linkS4class{NIFTIMetaInfo}}, \code{\linkS4class{AFNIMetaInfo}}
#'
#' @importFrom methods new
#' @importFrom assertthat assert_that
#'
#' @export
MetaInfo <- function(Dim, spacing, origin = rep(0, length(spacing)),
                    data_type = "FLOAT", label = "",
                    spatial_axes = OrientationList3D$AXIAL_LPI,
                    additional_axes = NullAxis) {

  assert_that(is.numeric(Dim) && all(Dim > 0) && all(Dim == floor(Dim)),
              msg = "'Dim' must be a vector of positive integers")

  assert_that(is.numeric(spacing) && all(spacing > 0),
              msg = "'spacing' must be a vector of positive numbers")

  assert_that(is.numeric(origin) && all(is.finite(origin)),
              msg = "'origin' must be a vector of finite numbers")

  # Validate data type
  valid_types <- c("BYTE", "SHORT", "INT", "FLOAT", "DOUBLE")
  assert_that(data_type %in% valid_types,
              msg = paste("'data_type' must be one of:", paste(valid_types, collapse = ", ")))

  # Create object
  new("MetaInfo",
      dims = as.integer(Dim),
      spacing = as.numeric(spacing),
      origin = as.numeric(origin),
      data_type = data_type,
      label = as.character(label),
      spatial_axes = spatial_axes,
      additional_axes = additional_axes)
}

#' Create NIFTIMetaInfo Object
#'
#' @title Create NIFTI Format Metadata Object
#' @description Creates a NIFTIMetaInfo object containing format-specific metadata
#'   for NIFTI format neuroimaging files.
#'
#' @param descriptor NIFTIFormat object specifying file format details
#' @param nifti_header List containing NIFTI header information
#'
#' @return A NIFTIMetaInfo object
#'
#' @details
#' The NIFTIMetaInfo object extends MetaInfo with NIFTI-specific features:
#' \itemize{
#'   \item NIFTI header fields (qform, sform matrices)
#'   \item Data scaling (slope, intercept)
#'   \item File organization (separate vs. single file)
#'   \item Orientation information
#' }
#'
#' Validation ensures:
#' \itemize{
#'   \item Valid NIFTI format
#'   \item Consistent dimensions
#'   \item Valid transformation matrices
#'   \item Proper data scaling
#' }
#'
#' @examples
#' \dontrun{
#' # Read NIFTI header
#' header <- readNIftiHeader("brain.nii")
#'
#' # Create format descriptor
#' fmt <- NIFTIFormat()
#'
#' # Create metadata
#' meta <- NIFTIMetaInfo(fmt, header)
#'
#' # Check dimensions
#' dim(meta)
#' }
#'
#' @seealso
#' \code{\link{MetaInfo}}
#'
#' @export
NIFTIMetaInfo <- function(descriptor, nifti_header) {
  # Validate inputs
  if (!inherits(descriptor, "NIFTIFormat")) {
    stop("'descriptor' must be a NIFTIFormat object")
  }
  if (!is.list(nifti_header)) {
    stop("'nifti_header' must be a list")
  }
  if (is.null(nifti_header$file_type) || tolower(nifti_header$file_type) != "nifti") {
    stop("Invalid NIFTI header: missing or incorrect file_type")
  }

  # Validate dimensions
  dims <- try(niftiDim(nifti_header), silent = TRUE)
  if (inherits(dims, "try-error")) {
    stop("Invalid dimensions in NIFTI header")
  }

  # Validate transformation
  if (!is.matrix(nifti_header$qform) || nrow(nifti_header$qform) != 4 ||
      ncol(nifti_header$qform) != 4) {
    stop("Invalid qform matrix in NIFTI header")
  }

  # Create object with validation
  tryCatch({
    new("NIFTIMetaInfo",
        header_file = header_file(descriptor, nifti_header$file_name),
        data_file = data_file(descriptor, nifti_header$file_name),
        descriptor = descriptor,
        endian = nifti_header$endian,
        data_offset = nifti_header$vox_offset,
        data_type = nifti_header$data_storage,
        bytes_per_element = as.integer(.getDataSize(nifti_header$data_storage)),
        dims = dims,
        spatial_axes = .nearestAnatomy(nifti_header$qform),
        additional_axes = NullAxis,
        spacing = nifti_header$pixdim[2:4],
        origin = nifti_header$qoffset,
        label = strip_extension(descriptor, basename(nifti_header$file_name)),
        intercept = nifti_header$scl_intercept,
        slope = nifti_header$scl_slope,
        header = nifti_header)
  }, error = function(e) {
    stop("Failed to create NIFTIMetaInfo: ", e$message)
  })
}

#' Create AFNIMetaInfo Object
#'
#' @title Create AFNI Format Metadata Object
#' @description Creates an AFNIMetaInfo object containing format-specific metadata
#'   for AFNI format neuroimaging files.
#'
#' @param descriptor AFNIFormat object specifying file format details
#' @param afni_header List containing AFNI header information
#'
#' @return An AFNIMetaInfo object
#'
#' @details
#' The AFNIMetaInfo object extends MetaInfo with AFNI-specific features:
#' \itemize{
#'   \item AFNI brick structure
#'   \item Sub-brick labels and scaling
#'   \item Space transformation
#'   \item Statistical parameters
#' }
#'
#' The function handles:
#' \itemize{
#'   \item Dimension extraction and validation
#'   \item Label generation for sub-bricks
#'   \item Transformation from AFNI to NIFTI space
#'   \item Data type and scaling setup
#' }
#'
#' @examples
#' \dontrun{
#' # Read AFNI header
#' header <- read_afni_header("brain+orig.HEAD")
#'
#' # Create format descriptor
#' fmt <- AFNIFormat()
#'
#' # Create metadata
#' meta <- AFNIMetaInfo(fmt, header)
#'
#' # Check dimensions
#' dim(meta)
#' }
#'
#'
#' @export
AFNIMetaInfo <- function(descriptor, afni_header) {
  # Validate inputs
  if (!inherits(descriptor, "AFNIFormat")) {
    stop("'descriptor' must be an AFNIFormat object")
  }
  if (!is.list(afni_header)) {
    stop("'afni_header' must be a list")
  }

  # Extract and validate dimensions
  if (is.null(afni_header$DATASET_DIMENSIONS) ||
      is.null(afni_header$DATASET_DIMENSIONS$content)) {
    stop("Missing DATASET_DIMENSIONS in AFNI header")
  }

  .Dim <- afni_header$DATASET_DIMENSIONS$content[
    afni_header$DATASET_DIMENSIONS$content > 0
  ]

  if (length(.Dim) < 3) {
    stop("AFNI dataset must have at least 3 dimensions")
  }

  # Add time dimension if present
  if (!is.null(afni_header$DATASET_RANK$content) &&
      afni_header$DATASET_RANK$content[2] > 1) {
    .Dim <- c(.Dim, afni_header$DATASET_RANK$content[2])
  }

  # Generate or extract labels
  labs <- if (is.null(afni_header$BRICK_LABS$content)) {
    paste0("#", seq(0, afni_header$DATASET_RANK$content[2] - 1))
  } else {
    afni_header$BRICK_LABS$content
  }

  # Calculate space transformation
  Tdicom <- tryCatch({
    matrix(afni_header$IJK_TO_DICOM$content, 3, 4, byrow = TRUE)
  }, error = function(e) {
    stop("Invalid IJK_TO_DICOM transformation in AFNI header")
  })

  TLPI <- perm_mat(OrientationList3D$AXIAL_RAI) %*% Tdicom[1:3, ]
  TLPI <- rbind(TLPI, c(0, 0, 0, 1))

  # Determine data type and size
  data_type <- switch(as.character(afni_header$BRICK_TYPES$content[1]),
                     "0" = "BYTE",
                     "1" = "SHORT",
                     "3" = "FLOAT",
                     stop("Unsupported BRICK_TYPE in AFNI header"))

  bytes_per_element <- switch(as.character(afni_header$BRICK_TYPES$content[1]),
                            "0" = 1L,
                            "1" = 2L,
                            "3" = 4L)

  # Create object with validation
  tryCatch({
    new("AFNIMetaInfo",
        header_file = header_file(descriptor, afni_header$file_name),
        data_file = data_file(descriptor, afni_header$file_name),
        descriptor = descriptor,
        endian = ifelse(afni_header[["BYTEORDER_STRING"]]$content == "MSB_FIRST",
                       "big", "little"),
        data_offset = 0L,
        data_type = data_type,
        bytes_per_element = as.integer(bytes_per_element),
        dims = .Dim,
        spatial_axes = .nearestAnatomy(TLPI),
        additional_axes = NullAxis,
        spacing = abs(afni_header$DELTA$content),
        origin = afni_header$ORIGIN$content,
        label = labs,
        intercept = 0,
        slope = ifelse(afni_header$BRICK_FLOAT_FACS$content == 0,
                      1,
                      afni_header$BRICK_FLOAT_FACS$content),
        header = afni_header)
  }, error = function(e) {
    stop("Failed to create AFNIMetaInfo: ", e$message)
  })
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

#' @rdname MetaInfo-methods
#' @aliases trans,MetaInfo-method
#'          trans,NIFTIMetaInfo-method
#'          data_reader,AFNIMetaInfo-method
#'          data_reader,NIFTIMetaInfo-method
#' @export
