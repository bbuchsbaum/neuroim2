% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/roi.R
\name{ROICoords}
\alias{ROICoords}
\title{Create ROI Coordinates Object}
\usage{
ROICoords(coords)
}
\arguments{
\item{coords}{A matrix with 3 columns representing (x, y, z) coordinates}
}
\value{
An \code{\linkS4class{ROICoords}} object
}
\description{
Creates an \code{\linkS4class{ROICoords}} object from a matrix of coordinates
representing points in 3D space.
}
\details{
ROI Coordinates
}
\examples{
coords <- matrix(c(1,2,3, 4,5,6), ncol=3, byrow=TRUE)
roi_coords <- ROICoords(coords)

}
