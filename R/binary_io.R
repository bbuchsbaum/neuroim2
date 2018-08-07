#' @include all_class.R
{}
#' @include all_generic.R
{}

.read_mmap <- function(meta, idx) {

  if (.Platform$endian != meta@endian) {
    stop(".read_mmap: swapped endian data not supported.")
    ## read raw bytes
    rawbytes <- mmap::mmap(meta@dataFile, mode=mmap::char(),
                           prot=mmap::mmapFlags("PROT_READ"))
    rawbytes <- rawbytes[(meta@dataOffset+1):length(rawbytes)]

    mmap::munmap(rawbytes)
    readBin(rawbytes, what=.getRStorage(meta@dataType), size=.getDataSize(meta@dataType), n=nels, endian=meta@endian)
  } else {
    #mmap::mmap(meta@dataFile, mode=.getMMapMode(meta@dataType), off=meta@dataOffset,prot=mmap::mmapFlags("PROT_READ"),
    #flags=mmap::mmapFlags("MAP_PRIVATE"))
    ret <- mmap::mmap(meta@data_file, mode=.getMMapMode(meta@data_type), prot=mmap::mmapFlags("PROT_READ"))
    offset <- meta@data_offset/.getDataSize(meta@data_type)
    idx_off <- idx + offset
    vals <- ret[idx_off]
    mmap::munmap(ret)
    vals
  }


}

read_mapped_series <- function(meta, idx) {
  if (endsWith(meta@data_file, ".gz")) {
    stop(paste("Cannot create series_reader with gzipped file", file_name))
  }

  assert_that(length(meta@dims) == 4, msg="'file_name' argument must refer to a 4-dimensional image")
  nels <- prod(meta@dims[1:3])

  dtype <- .getRStorage(meta@data_type)
  idx_set <- map(seq(1, meta@dims[4]), ~ idx + (nels*(.-1))) %>% flatten_dbl()
  ret <- .read_mmap(meta, idx_set)
  t(matrix(ret, length(idx), meta@dims[4]))
}

read_mapped_vols <- function(meta, idx) {
  if (endsWith(meta@data_file, ".gz")) {
    stop(paste("Cannot create series_reader with gzipped file", file_name))
  }

  assert_that(length(meta@dims) == 4, msg="'file_name' argument must refer to a 4-dimensional image")
  nels <- prod(meta@dims[1:3])
  nimages <- meta@dims[4]

  assert_that(min(idx) >= 1 && max(idx) <= nimages)

  dtype <- .getRStorage(meta@data_type)

  idx_set <- map(idx, ~ (.-1)*nels + seq(1,nels)) %>% flatten_dbl()
  ret <- .read_mmap(meta, idx_set)
  mat <- matrix(ret, nels, length(idx))
}


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
		new("BinaryReader", input=file(input, open="rb"), byte_offset=as.integer(byte_offset),
		    data_type=data_type, bytes_per_element=as.integer(bytes_per_element), endian=endian)
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
