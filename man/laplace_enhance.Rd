% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spat_filter.R
\name{laplace_enhance}
\alias{laplace_enhance}
\title{Laplacian Enhancement Filter for Volumetric Images}
\usage{
laplace_enhance(
  vol,
  mask,
  k = 2,
  patch_size = 3,
  search_radius = 2,
  h = 0.7,
  mapping_params = NULL,
  use_normalization_free = TRUE
)
}
\arguments{
\item{vol}{A \code{\linkS4class{NeuroVol}} object representing the image volume to be enhanced.}

\item{mask}{A \code{\linkS4class{LogicalNeuroVol}} object specifying the region to process. If not provided,
the entire volume will be processed.}

\item{k}{An integer specifying the number of layers in the decomposition (default is 2).}

\item{patch_size}{An integer specifying the size of patches for non-local means. Must be odd (default is 3).}

\item{search_radius}{An integer specifying the radius of the search window (default is 2).}

\item{h}{A numeric value controlling the filtering strength. Higher values mean more smoothing (default is 0.7).}

\item{mapping_params}{An optional list of parameters for the enhancement mappings.}

\item{use_normalization_free}{Logical indicating whether to use normalization-free weights (default is TRUE).}
}
\value{
A \code{\linkS4class{NeuroVol}} object representing the enhanced image.
}
\description{
This function applies a multi-layer Laplacian enhancement filter to a volumetric image (3D brain MRI data).
The filter enhances details while preserving edges using a non-local means approach with multiple scales.
}
