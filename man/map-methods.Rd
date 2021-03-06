% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R, R/neurovol.R
\name{mapf}
\alias{mapf}
\alias{mapf,NeuroVol,Kernel-method}
\title{apply a function to an object}
\usage{
mapf(x, m, ...)

\S4method{mapf}{NeuroVol,Kernel}(x, m, mask = NULL)
}
\arguments{
\item{x}{the object that is mapped}

\item{m}{the mapping object}

\item{...}{additional arguments}

\item{mask}{restrict application of kernel to masked area}
}
\description{
apply a function to an object

apply a kernel function to a \code{\linkS4class{NeuroVol}}
}
