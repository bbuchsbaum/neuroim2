% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/roi.R
\name{ROIVol}
\alias{ROIVol}
\title{Create an instance of class \code{\linkS4class{ROIVol}}}
\usage{
ROIVol(vspace, coords, data = rep(nrow(coords), 1))
}
\arguments{
\item{coords}{matrix of voxel coordinates}

\item{data}{the data values, numeric vector}

\item{space}{an instance of class \code{NeuroSpace}}
}
\value{
an instance of class \code{ROIVol}
}
\description{
Create an instance of class \code{\linkS4class{ROIVol}}
}