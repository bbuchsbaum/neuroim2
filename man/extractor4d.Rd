% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/array_like.R
\name{extractor4d}
\alias{extractor4d}
\alias{[,ArrayLike4D,matrix,missing,ANY-method}
\alias{[,ArrayLike4D,numeric,numeric,ANY-method}
\alias{[,ArrayLike4D,numeric,missing,ANY-method}
\alias{[,ArrayLike4D,missing,missing,ANY-method}
\alias{[,ArrayLike4D,missing,numeric,ANY-method}
\title{Array-like access for 4-dimensional data structures}
\usage{
\S4method{[}{ArrayLike4D,matrix,missing,ANY}(x, i, j, k, m, ..., drop = TRUE)

\S4method{[}{ArrayLike4D,numeric,numeric,ANY}(x, i, j, k, m, ..., drop = TRUE)

\S4method{[}{ArrayLike4D,numeric,missing,ANY}(x, i, j, k, m, ..., drop = TRUE)

\S4method{[}{ArrayLike4D,missing,missing,ANY}(x, i, j, k, m, ..., drop = TRUE)

\S4method{[}{ArrayLike4D,missing,numeric,ANY}(x, i, j, k, m, ..., drop = TRUE)
}
\arguments{
\item{x}{the object}

\item{i}{first index}

\item{j}{second index}

\item{k}{third index}

\item{m}{the fourth index}

\item{...}{additional args}

\item{drop}{dimension}
}
\description{
Array-like access for 4-dimensional data structures
}
