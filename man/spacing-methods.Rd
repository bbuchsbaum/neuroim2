% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R, R/neuro_obj.R, R/neurospace.R
\name{spacing}
\alias{spacing}
\alias{spacing,ROICoords-method}
\alias{spacing,NeuroObj-method}
\alias{spacing,NeuroSpace-method}
\title{extract the voxel dimensions of an image}
\usage{
spacing(x)

\S4method{spacing}{ROICoords}(x)

\S4method{spacing}{NeuroObj}(x)

\S4method{spacing}{NeuroSpace}(x)
}
\arguments{
\item{x}{the object}
}
\value{
a numeric vector
}
\description{
extract the voxel dimensions of an image

spacing
}
\examples{
bspace <- NeuroSpace(c(10,10,10), c(2,2,2))
all.equal(spacing(bspace), c(2,2,2))
}
