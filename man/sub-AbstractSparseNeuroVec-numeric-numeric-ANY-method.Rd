% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sparse_neurovec.R
\name{[,AbstractSparseNeuroVec,numeric,numeric,ANY-method}
\alias{[,AbstractSparseNeuroVec,numeric,numeric,ANY-method}
\title{Extractor Method for AbstractSparseNeuroVec}
\usage{
\S4method{[}{AbstractSparseNeuroVec,numeric,numeric,ANY}(x, i, j, k, m, ..., drop = TRUE)
}
\arguments{
\item{x}{An object of class \code{AbstractSparseNeuroVec}}

\item{i}{Numeric vector specifying the indices for the first dimension}

\item{j}{Numeric vector specifying the indices for the second dimension}

\item{k}{Numeric vector specifying the indices for the third dimension (optional)}

\item{m}{Numeric vector specifying the indices for the fourth dimension (optional)}

\item{...}{Additional arguments passed to methods}

\item{drop}{Logical indicating whether to drop dimensions of length one (default: TRUE)}
}
\value{
An array containing the extracted subset
}
\description{
Extracts a subset of data from a sparse four-dimensional brain image based on provided indices.
}
