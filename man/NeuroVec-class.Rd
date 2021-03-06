% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_class.R, R/neurovec.R
\docType{class}
\name{NeuroVec-class}
\alias{NeuroVec-class}
\alias{NeuroVec}
\title{NeuroVec}
\usage{
NeuroVec(data, space = NULL, mask = NULL, label = "")
}
\arguments{
\item{data}{the image data which can be a \code{matrix}, a 4d \code{array}, or a list of \code{NeuroVols}.
If the latter, the geometric space of the data \code{NeuroSpace} will be inferred from the constituent volumes,
which must all be identical.}

\item{space}{a \code{\linkS4class{NeuroSpace}} object. Does not ned to be included if \code{data} argument is a list of \code{NeuroVols}}

\item{mask}{an optional \code{array} of type \code{logical}}

\item{label}{a label of type \code{character}}
}
\value{
a concrete instance of \code{\linkS4class{NeuroVec}} class.
If \code{mask} is provided then \code{\linkS4class{SparseNeuroVec}}, otherwise \code{\linkS4class{DenseNeuroVec}}
}
\description{
Four-dimensional brain image

constructor function for virtual class \code{\linkS4class{NeuroVec}}
}
