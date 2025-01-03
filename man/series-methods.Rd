% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R, R/neurovec.R, R/neurovecseq.R,
%   R/sparse_neurovec.R
\name{series_roi}
\alias{series_roi}
\alias{series}
\alias{series,NeuroVec,matrix-method}
\alias{series_roi,NeuroVec,matrix-method}
\alias{series,NeuroVec,ROICoords-method}
\alias{series_roi,NeuroVec,ROICoords-method}
\alias{series,NeuroVec,LogicalNeuroVol-method}
\alias{series,NeuroVec,NeuroVol-method}
\alias{series_roi,NeuroVec,LogicalNeuroVol-method}
\alias{series,NeuroVec,integer-method}
\alias{series,DenseNeuroVec,integer-method}
\alias{series,NeuroVec,numeric-method}
\alias{series_roi,NeuroVec,numeric-method}
\alias{series,NeuroVecSeq,numeric-method}
\alias{series,NeuroVecSeq,matrix-method}
\alias{series_roi,NeuroVecSeq,matrix-method}
\alias{series,AbstractSparseNeuroVec,ROICoords-method}
\alias{series,AbstractSparseNeuroVec,matrix-method}
\alias{series,AbstractSparseNeuroVec,numeric-method}
\alias{series,AbstractSparseNeuroVec,integer-method}
\title{Extract one or more series from object and return as ROI object}
\usage{
series_roi(x, i, ...)

series(x, i, ...)

\S4method{series}{NeuroVec,matrix}(x, i)

\S4method{series_roi}{NeuroVec,matrix}(x, i)

\S4method{series}{NeuroVec,ROICoords}(x, i)

\S4method{series_roi}{NeuroVec,ROICoords}(x, i)

\S4method{series}{NeuroVec,LogicalNeuroVol}(x, i)

\S4method{series}{NeuroVec,NeuroVol}(x, i)

\S4method{series_roi}{NeuroVec,LogicalNeuroVol}(x, i)

\S4method{series}{NeuroVec,integer}(x, i, j, k, drop = TRUE)

\S4method{series}{DenseNeuroVec,integer}(x, i, j, k, drop = TRUE)

\S4method{series}{NeuroVec,numeric}(x, i, j, k)

\S4method{series_roi}{NeuroVec,numeric}(x, i, j, k)

\S4method{series_roi}{NeuroVec,numeric}(x, i, j, k)

\S4method{series_roi}{NeuroVec,matrix}(x, i)

\S4method{series_roi}{NeuroVec,LogicalNeuroVol}(x, i)

\S4method{series}{NeuroVecSeq,numeric}(x, i, j, k, drop = TRUE)

\S4method{series}{NeuroVecSeq,matrix}(x, i)

\S4method{series_roi}{NeuroVecSeq,matrix}(x, i)

\S4method{series}{AbstractSparseNeuroVec,ROICoords}(x, i)

\S4method{series}{AbstractSparseNeuroVec,matrix}(x, i)

\S4method{series}{AbstractSparseNeuroVec,numeric}(x, i, j, k)

\S4method{series}{AbstractSparseNeuroVec,integer}(x, i, j, k)
}
\arguments{
\item{x}{A NeuroVecSeq object}

\item{i}{A matrix of ROI coordinates (n x 3)}

\item{...}{additional arguments}

\item{j}{index for 2nd dimension}

\item{k}{index for 3rd dimension}

\item{drop}{whether to drop dimension of length 1}
}
\value{
A matrix where each column represents a voxel's time series

A matrix where each column represents a voxel's time series

A ROIVec object containing the time series for the specified ROI
}
\description{
Extract one or more series from object and return as ROI object

Extract one or more series from object

setMethod("series", signature(x="NeuroVecSeq", i="integer"),
          def=function(x, i, j, k, drop=TRUE) {
            if (missing(j) && missing(k)) {
              # Linear spatial indexing case
              # Get series from each vector and combine them
              series_list <- purrr::map(x@vecs, ~ series(., i))
              # Stack the matrices vertically
              do.call(rbind, series_list)
            } else {
              # Voxel coordinate indexing case
              assertthat::assert_that(length(i) == 1 && length(j) == 1 && length(k) == 1)
              purrr::map(x@vecs, ~ series(., i, j, k)) %>%
                purrr::flatten_dbl()
            }
          })
}
