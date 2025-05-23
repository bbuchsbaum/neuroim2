% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mapped_neurovec.R, R/neurovec.R, R/roi.R,
%   R/sparse_neurovec.R
\name{as.matrix,MappedNeuroVec-method}
\alias{as.matrix,MappedNeuroVec-method}
\alias{as.matrix,NeuroVec-method}
\alias{as.matrix,DenseNeuroVec-method}
\alias{as.matrix,ROIVec-method}
\alias{as.matrix,SparseNeuroVec-method}
\title{convert a \code{NeuroVec} to a matrix}
\usage{
\S4method{as.matrix}{MappedNeuroVec}(x)

\S4method{as.matrix}{NeuroVec}(x)

\S4method{as.matrix}{NeuroVec}(x)

\S4method{as.matrix}{DenseNeuroVec}(x)

\S4method{as.matrix}{ROIVec}(x)

\S4method{as.matrix}{SparseNeuroVec}(x, ...)
}
\arguments{
\item{x}{The object to convert to a matrix}

\item{...}{Additional arguments}
}
\value{
A matrix representation of the object
}
\description{
convert a \code{NeuroVec} to a matrix

convert a \code{\linkS4class{ROIVec}} to a matrix

Convert to Matrix
}
