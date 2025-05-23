% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R, R/common.R
\name{split_reduce}
\alias{split_reduce}
\alias{split_reduce,matrix,integer,function-method}
\alias{split_reduce,matrix,factor,missing-method}
\alias{split_reduce,matrix,factor,function-method}
\alias{split_reduce,NeuroVec,factor,function-method}
\alias{split_reduce,NeuroVec,factor,missing-method}
\title{Summarize Subsets of an Object by Splitting by Row and Applying a Summary Function}
\usage{
split_reduce(x, fac, FUN)

\S4method{split_reduce}{matrix,integer,function}(x, fac, FUN)

\S4method{split_reduce}{matrix,factor,missing}(x, fac)

\S4method{split_reduce}{matrix,factor,function}(x, fac, FUN)

\S4method{split_reduce}{NeuroVec,factor,function}(x, fac, FUN)

\S4method{split_reduce}{NeuroVec,factor,missing}(x, fac, FUN)
}
\arguments{
\item{x}{A numeric matrix or matrix-like object.}

\item{fac}{A factor to define subsets of the object.}

\item{FUN}{The summary function to apply to each subset. If not provided, the mean of each sub-matrix column is computed.}
}
\value{
A \code{matrix} (or matrix-like object) containing the summarized values after applying \code{FUN}.
}
\description{
This function summarizes subsets of a numeric matrix or matrix-like object by first splitting the object by row and then applying a summary function.
}
\details{
If 'FUN' is supplied, it must take a vector and return a single scalar value. If it returns more than one value, an error will occur.

If 'x' is a NeuroVec instance, voxels (dimensions 1:3) are treated as columns and time-series (dimension 4) as rows.
The summary function is then applied to groups of voxels. However, if the goal is to apply a function to groups of time-points.
}
\examples{
mat = matrix(rnorm(100*100), 100, 100)
fac = factor(sample(1:3, nrow(mat), replace=TRUE))
## Compute column means of each sub-matrix
ms <- split_reduce(mat, fac)
all.equal(row.names(ms), levels(fac))

## Compute column medians of each sub-matrix
ms <- split_reduce(mat, fac, median)

## Compute time-series means grouped over voxels.
## Here, 'length(fac)' must equal the number of voxels: 'prod(dim(bvec)[1:3])'
bvec <- NeuroVec(array(rnorm(24*24*24*24), c(24,24,24,24)), NeuroSpace(c(24,24,24,24), c(1,1,1)))
fac <- factor(sample(1:3, prod(dim(bvec)[1:3]), replace=TRUE))
ms <- split_reduce(bvec, fac)
ms2 <- split_reduce(bvec, fac, mean)
all.equal(row.names(ms), levels(fac))
all.equal(ms, ms2)

}
