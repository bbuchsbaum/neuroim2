% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R, R/neurospace.R
\name{dim_of}
\alias{dim_of}
\alias{dim_of,NeuroSpace,NamedAxis-method}
\title{Get the length of a given dimension of an object}
\usage{
dim_of(x, axis)

\S4method{dim_of}{NeuroSpace,NamedAxis}(x, axis)
}
\arguments{
\item{x}{The NeuroSpace object}

\item{axis}{The NamedAxis to query}
}
\value{
An integer representing the length of the specified \code{axis} of \code{x}.
}
\description{
This function returns the length of a given axis (dimension) of an object. The
axis can be specified using its position or name.
}
\examples{

x <- NeuroSpace(c(10,10,10), spacing=c(1,1,1))
stopifnot(dim_of(x, x@axes@i) == 10)
}
