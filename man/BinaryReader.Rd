% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binary_io.R
\name{BinaryReader}
\alias{BinaryReader}
\title{Create Binary Reader Object}
\usage{
BinaryReader(
  input,
  byte_offset,
  data_type,
  bytes_per_element,
  endian = .Platform$endian,
  signed = TRUE
)
}
\arguments{
\item{input}{Character string (file name) or connection object to read from}

\item{byte_offset}{Integer specifying bytes to skip at start of input}

\item{data_type}{Character string specifying R data type ('integer', 'double', etc.)}

\item{bytes_per_element}{Integer specifying bytes per data element (e.g., 4 or 8)}

\item{endian}{Character string specifying endianness ('big' or 'little', default: platform-specific)}

\item{signed}{Logical indicating if data type is signed (default: TRUE)}
}
\value{
An object of class \linkS4class{BinaryReader}
}
\description{
Create a new instance of the \linkS4class{BinaryReader} class for reading bulk binary data.
}
\examples{
\donttest{
# Create a temporary binary file
tmp <- tempfile()
writeBin(rnorm(100), tmp, size = 8)


# Read from existing connection with offset
con <- file(tmp, "rb")
reader <- BinaryReader(con, byte_offset=0,
                      data_type = "DOUBLE", bytes_per_element = 8L)
close(reader)

# Clean up
unlink(tmp)
}
}
\seealso{
\code{\link{BinaryWriter}} for writing binary data
}
