% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ops.R
\name{Compare,SparseNeuroVol,numeric-method}
\alias{Compare,SparseNeuroVol,numeric-method}
\alias{Compare,numeric,SparseNeuroVol-method}
\title{Compare a SparseNeuroVol object with a numeric value}
\usage{
\S4method{Compare}{SparseNeuroVol,numeric}(e1, e2)

\S4method{Compare}{numeric,SparseNeuroVol}(e1, e2)
}
\arguments{
\item{e1}{A numeric value to compare with the data of the SparseNeuroVol object.}

\item{e2}{A SparseNeuroVol object containing the data to be compared.}
}
\value{
The result of the comparison between the SparseNeuroVol object's data and the numeric value.

The result of the comparison between the numeric value and the SparseNeuroVol object's data.
}
\description{
This method compares the data of a SparseNeuroVol object (\code{e1}) with a numeric value (\code{e2}) using a generic comparison function.

This method compares a numeric value (\code{e1}) with the data of a SparseNeuroVol object (\code{e2}) using a generic comparison function.
}
