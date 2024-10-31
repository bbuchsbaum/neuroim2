#' @include all_class.R
{}
#' @include all_generic.R
{}

#' Read data using memory mapping
#'
#' This internal function reads data from a memory-mapped file based on the provided metadata and indices.
#'
#' @param meta An object containing metadata about the file to be read.
#' @param idx A vector of indices specifying which elements to read.
#' @return A vector of values read from the memory-mapped file.
#' @keywords internal
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

#' Read a mapped series from a 4D image file
#'
#' This function reads a series of data from a memory-mapped 4D image file.
#'
#' @param meta An object containing metadata about the file to be read.
#' @param idx A vector of indices specifying which elements to read.
#' @return A matrix of values read from the memory-mapped file, transposed.
#' @keywords internal
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

#' Read mapped data from a 4D image file
#'
#' This function reads data from a memory-mapped 4D image file.
#'
#' @param meta An object containing metadata about the file to be read.
#' @param idx A vector of indices specifying which elements to read.
#' @return A vector of values read from the memory-mapped file.
#' @keywords internal
read_mapped_data <- function(meta, idx) {
  if (endsWith(meta@data_file, ".gz")) {
    stop(paste("Cannot create series_reader with gzipped file", meta@data_file))
  }

  assert_that(length(meta@dims) == 4, msg="'file_name' argument must refer to a 4-dimensional image")
  nels <- prod(meta@dims[1:3])

  ret <- .read_mmap(meta, idx)
}

#' Read mapped volumes from a 4D image file
#'
#' This function reads volumes of data from a memory-mapped 4D image file.
#'
#' @param meta An object containing metadata about the file to be read.
#' @param idx A vector of indices specifying which volumes to read.
#' @return A matrix of values read from the memory-mapped file.
#' @keywords internal
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
}

#' Create a series reader function for a 4D image file
#'
#' This function creates a closure that can read series of data from a 4D image file.
#'
#' @param file_name The name of the 4D image file to read from.
#' @return A function that can read series of data from the specified file.
#' @keywords internal
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

#' Create a BinaryReader Object
#'
#' @description
#' This function creates a new instance of the \code{\linkS4class{BinaryReader}} class,
#' which supports reading of bulk binary data from a connection or file.
#'
#' @param input A file name (character string) to read from or a \code{connection} object.
#' @param byte_offset An integer specifying the number of bytes to skip at the start of input.
#' @param data_type A character string specifying the R data type of binary elements.
#' @param bytes_per_element An integer specifying the number of bytes in each data element 
#'   (e.g., 4 or 8 for floating point numbers).
#' @param endian A character string specifying the endianness of the binary input connection. 
#'   Default is \code{.Platform$endian}.
#' @param signed A logical value indicating whether the data type is signed (TRUE) or 
#'   unsigned (FALSE). Default is TRUE.
#'
#' @return A new instance of the \code{\linkS4class{BinaryReader}} class.
#'
#' @details
#' The \code{BinaryReader} function creates an object for efficient reading of binary data.
#' It can work with both file names and existing connection objects. If a file name is provided,
#' it opens a new connection in binary read mode.
#'
#' @seealso 
#' \code{\link{BinaryReader-class}} for the class definition and available methods.
#' \code{\link{BinaryWriter}} for writing binary data.
#'
#' @examples
#' \dontrun{
#' # Create a BinaryReader for a file
#' reader <- BinaryReader("data.bin", byte_offset = 0, data_type = "double", 
#'                        bytes_per_element = 8)
#'
#' # Create a BinaryReader for an existing connection
#' con <- file("data.bin", "rb")
#' reader <- BinaryReader(con, byte_offset = 100, data_type = "integer", 
#'                        bytes_per_element = 4, endian = "little")
#' }
#'
#' @export
#' @rdname BinaryReader
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

#' Create a BinaryWriter Object
#'
#' @description
#' This function creates a new instance of the \code{\linkS4class{BinaryWriter}} class,
#' which supports writing of bulk binary data to a connection or file.
#'
#' @param output A file name (character string) to write to or a \code{connection} object.
#' @param byte_offset An integer specifying the number of bytes to skip at the start of output.
#' @param data_type A character string specifying the R data type of binary elements to be written.
#' @param bytes_per_element An integer specifying the number of bytes in each data element 
#'   (e.g., 4 or 8 for floating point numbers).
#' @param endian A character string specifying the endianness of the binary output connection. 
#'   Default is \code{.Platform$endian}.
#'
#' @return A new instance of the \code{\linkS4class{BinaryWriter}} class.
#'
#' @details
#' The \code{BinaryWriter} function creates an object for efficient writing of binary data.
#' It can work with both file names and existing connection objects. If a file name is provided,
#' it opens a new connection in binary write mode.
#'
#' @seealso 
#' \code{\link{BinaryWriter-class}} for the class definition and available methods.
#' \code{\link{BinaryReader}} for reading binary data.
#'
#' @examples
#' \dontrun{
#' # Create a BinaryWriter for a file
#' writer <- BinaryWriter("output.bin", byte_offset = 0, data_type = "double", 
#'                        bytes_per_element = 8)
#'
#' # Create a BinaryWriter for an existing connection
#' con <- file("output.bin", "wb")
#' writer <- BinaryWriter(con, byte_offset = 100, data_type = "integer", 
#'                        bytes_per_element = 4, endian = "little")
#' }
#'
#' @export
#' @rdname BinaryWriter
BinaryWriter <- function(output, byte_offset, data_type, bytes_per_element, endian = .Platform$endian) {
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
setMethod(f="initialize", signature=signature("BinaryReader"),
		def=function(.Object, input, byte_offset, data_type, bytes_per_element, endian, signed) {
			.Object@input <- input
			.Object@byte_offset <- byte_offset
			.Object@data_type <- data_type
			.Object@bytes_per_element <- bytes_per_element
			.Object@endian <- endian
			.Object@signed <- signed

			## must be seekable connection, should enforce this
			##

			if (attr(.Object@input, "class")[[1]] != "gzfile") {
				seek(.Object@input, where=.Object@byte_offset, origin="start")
			} else {
				n <- as.integer(.Object@byte_offset/.Object@bytes_per_element)
				readBin(.Object@input, what=.Object@data_type, size=.Object@bytes_per_element, endian=.Object@endian, n=n)
			}

			.Object

		})

## code duplication, fix me
setMethod(f="initialize", signature=signature("BinaryWriter"),
		def=function(.Object, output, byte_offset, data_type, bytes_per_element, endian) {
			.Object@output <- output
			.Object@byte_offset <- byte_offset
			.Object@data_type <- data_type
			.Object@bytes_per_element <- bytes_per_element
			.Object@endian <- endian

			## must be seekable connection, should enforce this
			##
			#seek(.Object@output, where=.Object@byte_offset, origin="start")
			.Object
		})

#' Read Elements from a BinaryReader
#'
#' @description
#' This method reads a specified number of elements from a BinaryReader object.
#'
#' @param x A BinaryReader object.
#' @param num_elements A numeric value specifying the number of elements to read.
#' @return A vector of read elements.
#'
#' @export
#' @rdname read_elements-methods
setMethod(f="read_elements", signature=signature(x= "BinaryReader", num_elements="numeric"),
		def=function(x, num_elements) {
			readBin(x@input, what=x@data_type, size=x@bytes_per_element, n=num_elements, endian=x@endian, signed=x@signed)
		})

#' Write Elements to a BinaryWriter
#'
#' @description
#' This method writes elements to a BinaryWriter object.
#'
#' @param x A BinaryWriter object.
#' @param els A numeric vector of elements to write.
#'
#' @rdname write_elements-methods
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
#' @description
#' This method closes the connection associated with a BinaryReader or BinaryWriter object.
#'
#' @param con The BinaryReader or BinaryWriter object to close.
#'
#' @export
#' @rdname close-methods
setMethod(f="close", signature=signature(con= "BinaryReader"),
		def=function(con) {
			base::close(con@input)
		})


#' @rdname close-methods
setMethod(f="close", signature=signature(con= "BinaryWriter"),
		def=function(con) {
			base::close(con@output)
		})


#' Create a ColumnReader Object
#'
#' @description
#' This function creates a new instance of the \code{\linkS4class{ColumnReader}} class,
#' which represents a reader for column-oriented data.
#'
#' @param nrow An integer specifying the number of rows in the data.
#' @param ncol An integer specifying the number of columns in the data.
#' @param reader A function that takes a set of column indices and returns a matrix.
#'
#' @return A new instance of the \code{\linkS4class{ColumnReader}} class.
#'
#' @details
#' The \code{ColumnReader} function creates an object for efficient reading of column-oriented data.
#' The provided reader function should be capable of reading specified columns from the data source.
#'
#' @seealso 
#' \code{\link{ColumnReader-class}} for the class definition and available methods.
#'
#' @examples
#' \dontrun{
#' # Create a simple ColumnReader
#' reader_func <- function(cols) {
#'   # This is a placeholder function
#'   matrix(rnorm(100 * length(cols)), nrow = 100, ncol = length(cols))
#' }
#' col_reader <- ColumnReader(nrow = 100, ncol = 10, reader = reader_func)
#' }
#'
#' @export
#' @rdname ColumnReader
ColumnReader <- function(nrow, ncol, reader) {
  stopifnot(is.function(reader))
  new("ColumnReader", nrow=as.integer(nrow), ncol=as.integer(ncol), reader=reader)
}
