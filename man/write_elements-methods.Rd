% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R, R/binary_io.R
\name{write_elements}
\alias{write_elements}
\alias{write_elements,BinaryWriter,numeric-method}
\title{Write a sequence of elements from an input source}
\usage{
write_elements(x, els)

\S4method{write_elements}{BinaryWriter,numeric}(x, els)
}
\arguments{
\item{x}{A BinaryWriter object.}

\item{els}{A numeric vector of elements to write.}
}
\description{
This method writes elements to a BinaryWriter object.
}
