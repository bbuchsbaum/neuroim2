% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/neuroslice.R
\name{show,NeuroSlice-method}
\alias{show,NeuroSlice-method}
\title{Print NeuroSlice Object Details}
\usage{
\S4method{show}{NeuroSlice}(object)
}
\arguments{
\item{object}{A \code{NeuroSlice} object}
}
\description{
Displays a formatted summary of a \code{NeuroSlice} object's properties using
color-coded output for improved readability.
}
\details{
Display NeuroSlice Information


Shows key properties including:
\itemize{
  \item Object type and class
  \item Dimensions and memory usage
  \item Spatial properties (spacing, origin)
  \item Value range and statistics
  \item Axis orientations
}
}
\examples{
slice <- NeuroSlice(matrix(1:100, 10, 10), NeuroSpace(c(10, 10)))
show(slice)

}
