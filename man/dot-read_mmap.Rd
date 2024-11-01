% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binary_io.R
\name{.read_mmap}
\alias{.read_mmap}
\title{Read data using memory mapping}
\usage{
.read_mmap(meta, idx)
}
\arguments{
\item{meta}{An object containing metadata about the file to be read.}

\item{idx}{A vector of indices specifying which elements to read.}
}
\value{
A vector of values read from the memory-mapped file.
}
\description{
This internal function reads data from a memory-mapped file based on the provided metadata and indices.
}
\keyword{internal}
