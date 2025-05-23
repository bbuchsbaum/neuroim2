% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/neurovec.R
\name{NeuroVecSource}
\alias{NeuroVecSource}
\title{NeuroVecSource}
\usage{
NeuroVecSource(file_name, indices = NULL, mask = NULL)
}
\arguments{
\item{file_name}{The name of the 4-dimensional image file.}

\item{indices}{An optional integer vector specifying the subset of volume indices to load.
If not provided, all volumes will be loaded.}

\item{mask}{An optional logical array or \code{\linkS4class{NeuroVol}} object defining
the subset of voxels to load. If provided, a SparseNeuroVecSource object will be created.}
}
\value{
An instance of the \code{\linkS4class{NeuroVecSource}} class.
}
\description{
This function constructs a NeuroVecSource object, which represents the source
of a four-dimensional brain image.
}
\details{
If a \code{mask} is supplied, it should be a \code{\linkS4class{LogicalNeuroVol}} or
\code{\linkS4class{NeuroVol}} instance. If the latter, then the mask will be defined by
nonzero elements of the volume.
}
