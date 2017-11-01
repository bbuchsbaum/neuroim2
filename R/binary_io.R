#' @include all_class.R
NULL
#' @include all_generic.R
NULL



#' BinaryReader
#'
#' Constructor for  \code{\linkS4class{BinaryReader}} class
#'
#' @param input file name to read from or else a \code{connection} object
#' @param byte_offset the number of bytes to skip at the start of input
#' @param data_type R data type of binary elements
#' @param bytes_per_element number of bytes in each data element (e.g. 4 or 8 for floating point numbers)
#' @param endian endianness of binary input connection
#' @rdname BinaryReader
#' @export
BinaryReader <- function(input, byte_offset, data_type, bytes_per_element, endian=.Platform$endian) {
	if (is.character(input)) {
		new("BinaryReader", input=file(input, open="rb"), byte_offset=as.integer(byte_offset), data_type=data_type, bytes_per_element=as.integer(bytes_per_element), endian=endian)
	} else {
		stopifnot(inherits(input, "connection"))
		new("BinaryReader", input=input, byte_offset=as.integer(byte_offset), data_type=data_type, bytes_per_element=as.integer(bytes_per_element), endian=endian)
	}

}

#' BinaryWriter
#'
#' Constructor for  \code{\linkS4class{BinaryWriter}} class
#'
#' @param output file name to write to or else a \code{connection} object
#' @param byte_offset the number of bytes to skip at the start of output
#' @param data_type R data type of binary elements
#' @param bytes_per_element number of bytes in each data element (e.g. 4 or 8 for floating point numbers)
#' @param endian endianness of binary output connection
#' @rdname BinaryWriter-class
#' @export
BinaryWriter <- function(output, byte_offset, data_type, bytes_per_element, endian=.Platform$endian) {
	if (is.character(output)) {
		new("BinaryWriter", output=file(output, open="wb"), byte_offset=as.integer(byte_offset), data_type=data_type, bytes_per_element=as.integer(bytes_per_element), endian=endian)
	} else {
		stopifnot(inherits(output, "connection"))
		new("BinaryWriter", output=output, byte_offset=as.integer(byte_offset), data_type=data_type, bytes_per_element=as.integer(bytes_per_element), endian=endian)
	}

}

## code duplication, fix me. introduce "BinaryConnection superclass
setMethod(f="initialize", signature=signature("BinaryReader"),
		def=function(.Object, input, byte_offset, data_type, bytes_per_element, endian) {
			.Object@input <- input
			.Object@byte_offset <- byte_offset
			.Object@data_type <- data_type
			.Object@bytes_per_element <- bytes_per_element
			.Object@endian <- endian

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

#' read_elements
#'
#' @export
#' @rdname read_elements-methods
setMethod(f="read_elements", signature=signature(x= "BinaryReader", num_elements="numeric"),
		def=function(x, num_elements) {
			readBin(x@input, what=x@data_type, size=x@bytes_per_element, n=num_elements, endian=x@endian)
		})

#' write_elements
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

#' close
#'
#' @rdname close-methods
#' @export
#' @param con the object to close
setMethod(f="close", signature=signature(con= "BinaryReader"),
		def=function(con) {
			base::close(con@input)
		})


#' @rdname close-methods
setMethod(f="close", signature=signature(con= "BinaryWriter"),
		def=function(con) {
			base::close(con@output)
		})
