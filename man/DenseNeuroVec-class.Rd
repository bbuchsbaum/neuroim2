% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_class.R, R/neurovec.R
\docType{class}
\name{DenseNeuroVec-class}
\alias{DenseNeuroVec-class}
\alias{DenseNeuroVec}
\title{DenseNeuroVec}
\usage{
DenseNeuroVec(data, space, label = "")
}
\arguments{
\item{data}{a 4-dimensional \code{array} or a 2-dimension \code{matrix} that is either nvoxels by ntime-points or ntime-points by nvoxels}

\item{space}{a \code{\linkS4class{NeuroSpace}} object}

\item{label}{a label of type \code{character}}
}
\value{
\code{\linkS4class{DenseNeuroVec}} instance
}
\description{
Four-dimensional brain image, backed by an array

constructor function for class \code{\linkS4class{DenseNeuroVec}}
}
