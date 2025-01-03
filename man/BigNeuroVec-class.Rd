% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_class.R
\docType{class}
\name{BigNeuroVec-class}
\alias{BigNeuroVec-class}
\title{BigNeuroVec Class}
\description{
A class representing a sparse four-dimensional brain image backed by a disk-based
big matrix. BigNeuroVec objects are designed for efficient handling of large-scale
brain imaging data that exceeds available memory.
}
\details{
BigNeuroVec leverages file-backed storage to manage large 4D neuroimaging datasets
that would typically exceed available RAM. It combines the sparse representation
framework of \code{\linkS4class{AbstractSparseNeuroVec}} with the disk-based
storage capabilities of \code{FBM}, allowing for out-of-core computations on
massive datasets.
}
\section{Slots}{

\describe{
\item{\code{data}}{An instance of class \code{FBM} from the \code{bigstatsr} package,
containing time-series data. The FBM (File-Backed Big Matrix) is a matrix-like
structure stored on disk, enabling efficient handling of large-scale data.}
}}

\section{Inheritance}{

\code{BigNeuroVec} inherits from:
\itemize{
  \item \code{\linkS4class{NeuroVec}}: Base class for 4D brain images
  \item \code{\linkS4class{AbstractSparseNeuroVec}}: Provides sparse representation framework
  \item \code{ArrayLike4D}: Interface for 4D array-like operations
}
}

\examples{
\dontrun{
# Create a BigNeuroVec object
library(bigstatsr)

# Create a file-backed big matrix
fbm <- FBM(10000, 1000, init = rnorm(10000000))

# Create a mask for sparse representation
mask <- LogicalNeuroVol(array(runif(100*100*100) > 0.7, dim=c(100,100,100)))

# Create a BigNeuroVec object
big_vec <- BigNeuroVec(data = fbm, mask = mask, space = NeuroSpace(dim=c(100,100,100)))

# Access a subset of the data
subset <- big_vec[,,, 1:10]
}

}
\seealso{
\code{\link{AbstractSparseNeuroVec-class}} for the parent sparse representation class.
\code{\link{NeuroVec-class}} for the base 4D brain image class.
\code{\link[bigstatsr]{FBM}} for details on File-Backed Big Matrix objects.
}
