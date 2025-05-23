% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R, R/meta_info.R, R/neuro_obj.R,
%   R/neurospace.R
\name{trans}
\alias{trans}
\alias{trans,MetaInfo-method}
\alias{trans,NIFTIMetaInfo-method}
\alias{trans,NeuroObj-method}
\alias{trans,NeuroSpace-method}
\title{Extract image coordinate transformation}
\usage{
trans(x)

\S4method{trans}{MetaInfo}(x)

\S4method{trans}{NeuroObj}(x)

\S4method{trans}{NeuroSpace}(x)
}
\arguments{
\item{x}{an object with a transformation}
}
\value{
A numeric 4x4 matrix that maps from grid coordinates to real-world coordinates.
}
\description{
Extract image coordinate transformation

Get transformation matrix
}
\details{
This function returns a transformation that can be used to go from "grid coordinates" to "real world coordinates" in millimeters.
see \code{\linkS4class{NeuroSpace}}
}
\examples{
bspace <- NeuroSpace(c(10,10,10), c(2,2,2))
trans(bspace)
all.equal(dim(trans(bspace)), c(4,4))
}
