% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_class.R, R/binary_io.R
\docType{class}
\name{BinaryWriter-class}
\alias{BinaryWriter-class}
\alias{BinaryWriter}
\title{BinaryWriter}
\usage{
BinaryWriter(
  output,
  byte_offset,
  data_type,
  bytes_per_element,
  endian = .Platform$endian
)
}
\arguments{
\item{output}{file name to write to or else a \code{connection} object}

\item{byte_offset}{the number of bytes to skip at the start of output}

\item{data_type}{R data type of binary elements}

\item{bytes_per_element}{number of bytes in each data element (e.g. 4 or 8 for floating point numbers)}

\item{endian}{endianness of binary output connection}
}
\value{
a new instance of type \code{BinaryWriter}
}
\description{
This class supports writing of bulk binary data to a connection

Constructor for  \code{\linkS4class{BinaryWriter}} class
}
\section{Slots}{

\describe{
\item{\code{output}}{the binary output connection}

\item{\code{byte_offset}}{the number of bytes to skip at the start of input}

\item{\code{data_type}}{the dataType of the binary Elements}

\item{\code{bytes_per_element}}{number of bytes in each data element (e.g. 4 or 8 for floating point numbers)}

\item{\code{endian}}{endianness of binary output connection}
}}

