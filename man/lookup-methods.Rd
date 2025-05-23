% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/index_vol.R, R/sparse_neurovec.R
\name{lookup,IndexLookupVol,numeric-method}
\alias{lookup,IndexLookupVol,numeric-method}
\alias{lookup,AbstractSparseNeuroVec,numeric-method}
\title{Lookup Values in an IndexLookupVol Object}
\usage{
\S4method{lookup}{IndexLookupVol,numeric}(x, i)

\S4method{lookup}{AbstractSparseNeuroVec,numeric}(x, i)
}
\arguments{
\item{x}{An \code{\linkS4class{IndexLookupVol}} object}

\item{i}{A numeric vector of indices to look up}
}
\description{
Performs a lookup operation on an IndexLookupVol object.
}
\examples{

space <- NeuroSpace(c(64, 64, 64), c(1, 1, 1), c(0, 0, 0))
ilv <- IndexLookupVol(space, c(1:100))
lookup(ilv, c(1, 2, 3))  # Look up values for indices 1, 2, and 3


}
\keyword{internal}
