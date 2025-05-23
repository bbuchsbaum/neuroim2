% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R, R/filebacked_neurovec.R,
%   R/neurovec.R, R/neurovecseq.R, R/sparse_neurovec.R
\name{sub_vector}
\alias{sub_vector}
\alias{sub_vector,FileBackedNeuroVec,numeric-method}
\alias{sub_vector,NeuroVec,numeric-method}
\alias{sub_vector,NeuroVecSeq,numeric-method}
\alias{sub_vector,SparseNeuroVec,numeric-method}
\title{Generic function to extract a sub-vector from a \code{NeuroVec} object.}
\usage{
sub_vector(x, i, ...)

\S4method{sub_vector}{FileBackedNeuroVec,numeric}(x, i)

\S4method{sub_vector}{NeuroVec,numeric}(x, i)

\S4method{sub_vector}{NeuroVecSeq,numeric}(x, i)

\S4method{sub_vector}{NeuroVecSeq,numeric}(x, i)

\S4method{sub_vector}{SparseNeuroVec,numeric}(x, i)
}
\arguments{
\item{x}{A NeuroVecSeq object}

\item{i}{Numeric vector of indices specifying the time points to extract}

\item{...}{additional arguments}
}
\value{
A \code{NeuroVec} object that is a sub-sequence of the supplied object.

A NeuroVecSeq object containing the extracted subsequence
}
\description{
Extracts a subset of volumes from a file-backed neuroimaging vector and returns
them as a dense (in-memory) vector.

Extracts a subsequence of volumes from a NeuroVecSeq object.
}
\details{
This method efficiently reads only the requested volumes from disk, converting them
to an in-memory representation. The spatial metadata is preserved but adjusted to
reflect the new number of volumes.

Memory usage is proportional to the number of volumes requested, not the size of
the full dataset.
}
\examples{
bvec <- NeuroVec(array(rnorm(24*24*24*24), c(24,24,24,24)), NeuroSpace(c(24,24,24,24), c(1,1,1)))
vec <- sub_vector(bvec,1:2)
all.equal(2, dim(vec)[4])

vec <- sub_vector(bvec, c(1,3,5,7))
all.equal(4, dim(vec)[4])

mask <- LogicalNeuroVol(rep(TRUE, 24*24*24), NeuroSpace(c(24,24,24), c(1,1,1)))
svec <- SparseNeuroVec(array(rnorm(24*24*24*24), c(24,24,24,24)),
NeuroSpace(c(24,24,24,24), c(1,1,1)), mask)
vec <- sub_vector(svec, c(1,3,5))
all.equal(3, dim(vec)[4])
}
