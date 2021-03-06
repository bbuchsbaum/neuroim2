% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_class.R, R/neurovol.R
\docType{class}
\name{DenseNeuroVol-class}
\alias{DenseNeuroVol-class}
\alias{DenseNeuroVol}
\title{DenseNeuroVol}
\usage{
DenseNeuroVol(data, space, label = "", indices = NULL)
}
\arguments{
\item{data}{a three-dimensional \code{array}}

\item{space}{an instance of class \code{\linkS4class{NeuroSpace}}}

\item{label}{a \code{character} string}

\item{indices}{an optional 1-d index vector}
}
\value{
\code{\linkS4class{DenseNeuroVol}} instance
}
\description{
Three-dimensional brain image, backed by an \code{array}

Construct a \code{\linkS4class{DenseNeuroVol}} instance
}
