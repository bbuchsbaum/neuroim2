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
\item{x}{A FileFormat object.}

\item{file_name}{A character string representing the file name to be matched.}
}
\value{
TRUE for match, FALSE otherwise

A logical value: TRUE if the file name matches the header file format, 
  FALSE otherwise.
}
\description{
This method checks if a given file name conforms to the header file format
specified in the FileFormat object.
}
\details{
The function uses a regular expression to check if the file name ends with 
the header extension specified in the FileFormat object. It returns TRUE if 
there's a match, and FALSE otherwise.
}
\examples{
\dontrun{
file_format <- new("FileFormat", header_extension = "hdr", data_extension = "img")
header_file_matches(file_format, "brain_scan.hdr")
header_file_matches(file_format, "brain_scan.img")
}

}
\seealso{
\code{\link{FileFormat-class}}, \code{\link{file_matches}}, \code{\link{data_file_matches}}
}
