% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R, R/file_format.R
\name{header_file_matches}
\alias{header_file_matches}
\alias{header_file_matches,FileFormat,character-method}
\title{Generic function to test whether a file name conforms to the given \code{\linkS4class{FileFormat}} instance.
Will test for match to header file only}
\usage{
header_file_matches(x, file_name)

\S4method{header_file_matches}{FileFormat,character}(x, file_name)
}
\arguments{
\item{x}{A \linkS4class{FileFormat} object specifying the format requirements}

\item{file_name}{A character string specifying the file name to validate}
}
\value{
\code{TRUE} for match, \code{FALSE} otherwise.

A logical value: \code{TRUE} if the file name matches the header format,
  \code{FALSE} otherwise
}
\description{
Validates whether a file name conforms to the header file format specification.
}
\details{
The function performs case-sensitive pattern matching to verify that the file name
ends with the specified header extension. The match is performed using a regular
expression that ensures the extension appears at the end of the file name.
}
\examples{

fmt <- new("FileFormat", header_extension = "hdr", data_extension = "img")
header_file_matches(fmt, "brain_scan.hdr")  # TRUE
header_file_matches(fmt, "brain_scan.img")  # FALSE
header_file_matches(fmt, "brain.hdr.gz")    # FALSE


}
\seealso{
\code{\link{file_matches}}, \code{\link{data_file_matches}} for related
file format validation
}
