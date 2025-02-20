#' @include all_class.R
{}
#' @include all_generic.R
{}

## Helper Functions for enforcing seekability

#' Ensure that the input connection for BinaryReader is seekable
#'
#' @param con A connection object used for reading.
#' @param byte_offset Byte offset to move to.
#' @param bytes_per_element Number of bytes per element.
#' @param data_type The R data type.
#' @param endian The endianness of the data.
#'
#' @keywords internal
#' @noRd
ensure_reader_seekable <- function(con, byte_offset, bytes_per_element, data_type, endian) {
  if (!inherits(con, "gzfile")) {
    test <- try(seek(con, where = byte_offset, origin="start"), silent=TRUE)
    if(inherits(test, "try-error"))
      stop("Input connection must be seekable for Binary I/O.")
  } else {
    n <- as.integer(byte_offset / bytes_per_element)
    readBin(con, what = data_type, size = bytes_per_element, endian = endian, n = n)
  }
  invisible(NULL)
}

#' Ensure that the output connection for BinaryWriter is seekable
#'
#' @param con A connection object used for writing.
#' @param byte_offset Byte offset to move to.
#'
#' @keywords internal
#' @noRd
ensure_writer_seekable <- function(con, byte_offset) {
  if (!inherits(con, "gzfile")) {
    test <- try(seek(con, where = byte_offset, origin="start"), silent=TRUE)
    if(inherits(test, "try-error"))
      stop("Output connection must be seekable for Binary I/O.")
  } else {
    stop("Cannot use gzipped connection for BinaryWriter")
  }
  invisible(NULL)
}

#' Read Data Using Memory Mapping
#'
#' Read data from a memory-mapped file based on metadata and indices.
#'
#' @param meta An object of class \linkS4class{ImageMetadata} containing file metadata
#' @param idx Integer vector of indices specifying elements to read
#' @return A numeric vector of values read from the memory-mapped file
#' @keywords internal
#' @noRd
.read_mmap <- function(meta, idx) {
  if (.Platform$endian != meta@endian) {
    stop(".read_mmap: swapped endian data not supported.")
  }

  ret <- mmap::mmap(meta@data_file, mode=.getMMapMode(meta@data_type), prot=mmap::mmapFlags("PROT_READ"))
  offset <- meta@data_offset/.getDataSize(meta@data_type)
  idx_off <- idx + offset
  vals <- ret[idx_off]
  mmap::munmap(ret)
  vals
}

#' Read Mapped Series from 4D Image
#'
#' Read a time series of data from a memory-mapped 4D image file.
#'
#' @param meta An object of class \linkS4class{ImageMetadata} containing file metadata
#' @param idx Integer vector of indices specifying elements to read
#' @return A numeric matrix of values with dimensions [time, voxels]
#' @keywords internal
#' @noRd
read_mapped_series <- function(meta, idx) {
  if (endsWith(meta@data_file, ".gz")) {
    stop(paste("Cannot create series_reader with gzipped file", meta@data_file))
  }

  assert_that(length(meta@dims) == 4, msg="'file_name' argument must refer to a 4-dimensional image")
  nels <- prod(meta@dims[1:3])

  dtype <- .getRStorage(meta@data_type)
  idx_set <- map(seq(1, meta@dims[4]), ~ idx + (nels*(.-1))) %>% flatten_dbl()
  ret <- .read_mmap(meta, idx_set)
  t(matrix(ret, length(idx), meta@dims[4]))
}

#' Read Mapped Data from 4D Image
#'
#' Read data from a memory-mapped 4D image file.
#'
#' @param meta An object of class \linkS4class{ImageMetadata} containing file metadata
#' @param idx Integer vector of indices specifying elements to read
#' @return A numeric vector of values
#' @keywords internal
#' @noRd
read_mapped_data <- function(meta, idx) {
  if (endsWith(meta@data_file, ".gz")) {
    stop(paste("Cannot create series_reader with gzipped file", meta@data_file))
  }

  assert_that(length(meta@dims) == 4, msg="'file_name' argument must refer to a 4-dimensional image")
  nels <- prod(meta@dims[1:3])

  ret <- .read_mmap(meta, idx)
}

#' Read Mapped Volumes from 4D Image
#'
#' Read volume data from a memory-mapped 4D image file.
#'
#' @param meta An object of class \linkS4class{ImageMetadata} containing file metadata
#' @param idx Integer vector of indices specifying volumes to read
#' @return A numeric matrix with dimensions [time, voxels]
#' @keywords internal
#' @noRd
read_mapped_vols <- function(meta, idx) {
  if (endsWith(meta@data_file, ".gz")) {
    stop(paste("Cannot create series_reader with gzipped file", meta@data_file))
  }

  assert_that(length(meta@dims) == 4, msg="'file_name' argument must refer to a 4-dimensional image")
  nels <- prod(meta@dims[1:3])
  nimages <- meta@dims[4]

  assert_that(min(idx) >= 1 && max(idx) <= nimages)

  idx_set <- map(idx, ~ (.-1)*nels + seq(1,nels)) %>% flatten_dbl()
  ret <- .read_mmap(meta, idx_set)
  mat <- matrix(ret, nels, length(idx))
  t(mat)  # Transpose to get [time, voxels]
}

#' Create Series Reader for 4D Image
#'
#' Create a closure function that can read time series data from a 4D image file.
#'
#' @param file_name Character string specifying the path to the 4D image file
#' @return A function that accepts parameters:
#'   \itemize{
#'     \item first: Integer index of first element to read
#'     \item last: Integer index of last element to read
#'     \item input: Optional connection object (if NULL, creates new connection)
#'   }
#' @keywords internal
#' @noRd
series_reader <- function(file_name) {
  if (endsWith(file_name, ".gz")) {
    stop(paste("Cannot create series_reader with gzipped file", file_name))
  }

  meta <- read_header(file_name)
  assert_that(length(meta@dims) == 4, msg="'file_name' argument must refer to a 4-dimensional image")
  nels <- prod(meta@dims[1:3])

  dtype <- .getRStorage(meta@data_type)

  f <- function(first, last, input=NULL) {
    if (is.null(input)) {
      input <- file(file_name, open="rb")
      close <- TRUE
    } else {
      close <- FALSE
    }

    res <- lapply(seq(1, meta@dims[4]), function(i) {
      seek(input, where=(i-1)*nels*meta@bytes_per_element + meta@data_offset + (first-1)*meta@bytes_per_element, origin="start")
      readBin(input, what=dtype, size=meta@bytes_per_element, n=(last-first)+1, endian=meta@endian)
    })

    if (close) {
      close(input)
    }

    do.call(rbind, res)
  }

  f
}

#' Create Binary Reader Object
#'
#' Create a new instance of the \linkS4class{BinaryReader} class for reading bulk binary data.
#'
#' @param input Character string (file name) or connection object to read from
#' @param byte_offset Integer specifying bytes to skip at start of input
#' @param data_type Character string specifying R data type ('integer', 'double', etc.)
#' @param bytes_per_element Integer specifying bytes per data element (e.g., 4 or 8)
#' @param endian Character string specifying endianness ('big' or 'little', default: platform-specific)
#' @param signed Logical indicating if data type is signed (default: TRUE)
#' @return An object of class \linkS4class{BinaryReader}
#' @examples
#' \donttest{
#' # Create a temporary binary file
#' tmp <- tempfile()
#' writeBin(rnorm(100), tmp, size = 8)
#'
#'
#' # Read from existing connection with offset
#' con <- file(tmp, "rb")
#' reader <- BinaryReader(con, byte_offset=0,
#'                       data_type = "DOUBLE", bytes_per_element = 8L)
#' close(reader)
#'
#' # Clean up
#' unlink(tmp)
#' }
#' @seealso \code{\link{BinaryWriter}} for writing binary data
#' @export
BinaryReader <- function(input, byte_offset, data_type, bytes_per_element,
                        endian = .Platform$endian, signed = TRUE) {
  if (is.character(input)) {
    new("BinaryReader", input = file(input, open = "rb"),
        byte_offset = as.integer(byte_offset),
        data_type = data_type,
        bytes_per_element = as.integer(bytes_per_element),
        endian = endian, signed = signed)
  } else {
    stopifnot(inherits(input, "connection"))
    new("BinaryReader", input = input,
        byte_offset = as.integer(byte_offset),
        data_type = data_type,
        bytes_per_element = as.integer(bytes_per_element),
        endian = endian, signed = signed)
  }
}

#' Create Binary Writer Object
#'
#' Create a new instance of the \linkS4class{BinaryWriter} class for writing bulk binary data.
#'
#' @param output Character string (file name) or connection object to write to
#' @param byte_offset Integer specifying bytes to skip at start of output
#' @param data_type Character string specifying R data type ('integer', 'double', etc.)
#' @param bytes_per_element Integer specifying bytes per data element (e.g., 4 or 8)
#' @param endian Character string specifying endianness ('big' or 'little', default: platform-specific)
#' @return An object of class \linkS4class{BinaryWriter}
#' @examples
#' \donttest{
#'
#' tmp <- tempfile()
#' # Write to existing connection with offset
#' con <- file(tmp, "wb")
#' writer <- BinaryWriter(con, byte_offset = 100L,
#'                       data_type = "integer", bytes_per_element = 4L)
#' unlink(tmp)
#' }
#' @seealso \code{\link{BinaryReader}} for reading binary data
#' @export
BinaryWriter <- function(output, byte_offset, data_type, bytes_per_element,
                        endian = .Platform$endian) {
  if (is.character(output)) {
    new("BinaryWriter", output = file(output, open = "wb"),
        byte_offset = as.integer(byte_offset),
        data_type = data_type,
        bytes_per_element = as.integer(bytes_per_element),
        endian = endian)
  } else {
    stopifnot(inherits(output, "connection"))
    new("BinaryWriter", output = output,
        byte_offset = as.integer(byte_offset),
        data_type = data_type,
        bytes_per_element = as.integer(bytes_per_element),
        endian = endian)
  }
}

## some code duplication here.
setMethod(f="initialize", signature=signature(.Object="BinaryReader"),
		def=function(.Object, input, byte_offset, data_type, bytes_per_element, endian, signed) {
			.Object@input <- input
			.Object@byte_offset <- as.integer(byte_offset)
			.Object@data_type <- data_type
			.Object@bytes_per_element <- as.integer(bytes_per_element)
			.Object@endian <- endian
			.Object@signed <- signed

			# Ensure the reader's connection is seekable (or properly skipped if gzipped)
			ensure_reader_seekable(.Object@input, .Object@byte_offset, .Object@bytes_per_element, .Object@data_type, .Object@endian)

			.Object
		})

## code duplication, fix me
setMethod(f="initialize", signature=signature(.Object="BinaryWriter"),
		def=function(.Object, output, byte_offset, data_type, bytes_per_element, endian) {
			.Object@output <- output
			.Object@byte_offset <- as.integer(byte_offset)
			.Object@data_type <- data_type
			.Object@bytes_per_element <- as.integer(bytes_per_element)
			.Object@endian <- endian

			# Ensure that the writer's connection is seekable
			ensure_writer_seekable(.Object@output, .Object@byte_offset)

			.Object
		})

#' Read Elements from Binary Reader
#'
#' Read a specified number of elements from a \linkS4class{BinaryReader} object.
#'
#' @param x Object of class \linkS4class{BinaryReader}
#' @param num_elements Integer specifying number of elements to read
#' @return Numeric vector of read elements
#' @examples
#' \donttest{
#' # Create a temporary binary file with some test data
#' tmp <- tempfile()
#' con <- file(tmp, "wb")
#' test_data <- rnorm(100)
#' writeBin(test_data, con, size = 8)
#' close(con)
#'
#' # Create reader and read the data
#' reader <- BinaryReader(tmp, byte_offset = 0L,
#'                       data_type = "DOUBLE", bytes_per_element = 8L)
#' data <- read_elements(reader, 100)
#' close(reader)
#'
#' # Clean up
#' unlink(tmp)
#' }
#' @export
setMethod(f="read_elements", signature=signature(x= "BinaryReader", num_elements="numeric"),
		def=function(x, num_elements) {
			readBin(x@input, what=x@data_type, size=x@bytes_per_element, n=num_elements, endian=x@endian, signed=x@signed)
		})

#' Write Elements to Binary Writer
#'
#' Write elements to a \linkS4class{BinaryWriter} object.
#'
#' @param x Object of class \linkS4class{BinaryWriter}
#' @param els Numeric vector of elements to write
#' @examples
#' \donttest{
#' # Create a temporary binary file for writing
#' tmp <- tempfile()
#' writer <- BinaryWriter(tmp, byte_offset = 0L,
#'                       data_type = "DOUBLE", bytes_per_element = 8L)
#' # Write some data
#' write_elements(writer, rnorm(100))
#' close(writer)
#'
#' # Clean up
#' unlink(tmp)
#' }
#' @export
setMethod(f="write_elements", signature=signature(x= "BinaryWriter", els="numeric"),
		def=function(x, els) {
			if (.getRStorage(x@data_type) == "integer") {
				writeBin(as.integer(els), x@output, size=x@bytes_per_element, endian=x@endian)
			} else if (.getRStorage(x@data_type) == "double") {
				writeBin(as.double(els), x@output, size=x@bytes_per_element, endian=x@endian)
			} else {
				stop("failed to handle data type: ", x@data_type)
			}

		})

#' Close a BinaryReader or BinaryWriter
#'
#' Closes the underlying connection associated with a BinaryReader or BinaryWriter object.
#' This should be called when you're done with the reader/writer to free system resources.
#'
#' @param con The BinaryReader or BinaryWriter object to close.
#' @examples
#' \donttest{
#' # Create a temporary file and write some data
#' tmp <- tempfile()
#' writer <- BinaryWriter(tmp, byte_offset = 0L,
#'                       data_type = "DOUBLE", bytes_per_element = 8L)
#' write_elements(writer, rnorm(100))
#' close(writer)
#'
#' # Read the data back
#' reader <- BinaryReader(tmp, byte_offset = 0L,
#'                       data_type = "DOUBLE", bytes_per_element = 8L)
#' data <- read_elements(reader, 100)
#' close(reader)
#'
#' # Clean up
#' unlink(tmp)
#' }
#'
#' @export
#' @rdname close-methods
setMethod(f="close", signature=signature(con= "BinaryReader"),
		def=function(con) {
			base::close(con@input)
		})


#' @rdname close-methods
#' @export
setMethod(f="close", signature=signature(con= "BinaryWriter"),
		def=function(con) {
			base::close(con@output)
		})


#' Create Column Reader Object
#'
#' Create a new instance of the \linkS4class{ColumnReader} class for reading column-oriented data.
#'
#' @param nrow Integer specifying number of rows in data
#' @param ncol Integer specifying number of columns in data
#' @param reader Function that takes column indices and returns matrix
#' @return An object of class \linkS4class{ColumnReader}
#' @examples
#'
#' reader_func <- function(cols) {
#'   matrix(rnorm(100 * length(cols)), 100, length(cols))
#' }
#' col_reader <- ColumnReader(nrow = 100L, ncol = 10L, reader = reader_func)
#'
#' @export
ColumnReader <- function(nrow, ncol, reader) {
  stopifnot(is.function(reader))
  new("ColumnReader", nrow=as.integer(nrow), ncol=as.integer(ncol), reader=reader)
}
