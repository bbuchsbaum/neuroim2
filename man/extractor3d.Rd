% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/array_like.R, R/neurovol.R
\name{extractor3d}
\alias{extractor3d}
\alias{[,ArrayLike3D,numeric,missing,ANY-method}
\alias{[,ArrayLike3D,matrix,missing,ANY-method}
\alias{[,ArrayLike3D,missing,missing,ANY-method}
\alias{[,ArrayLike3D,missing,numeric,ANY-method}
\alias{[,NeuroVol,ROICoords,missing,ANY-method}
\title{Array-like access for 3-dimensional data structures}
\usage{
\S4method{[}{ArrayLike3D,numeric,missing,ANY}(x, i, j, k, ..., drop = TRUE)

\S4method{[}{ArrayLike3D,matrix,missing,ANY}(x, i, j, k, ..., drop = TRUE)

\S4method{[}{ArrayLike3D,missing,missing,ANY}(x, i, j, k, ..., drop = TRUE)

\S4method{[}{ArrayLike3D,missing,numeric,ANY}(x, i, j, k, ..., drop = TRUE)

\S4method{[}{NeuroVol,ROICoords,missing,ANY}(x, i, j, k, ..., drop = TRUE)
}
\arguments{
\item{x}{the object}

\item{i}{first index}

\item{j}{second index}

\item{k}{third index}

\item{...}{additional args}

\item{drop}{dimension}
}
\description{
Array-like access for 3-dimensional data structures
}
