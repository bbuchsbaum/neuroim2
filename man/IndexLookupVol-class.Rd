% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_class.R, R/index_vol.R
\docType{class}
\name{IndexLookupVol-class}
\alias{IndexLookupVol-class}
\alias{IndexLookupVol}
\title{IndexLookupVol Class}
\usage{
IndexLookupVol(space, indices)
}
\arguments{
\item{space}{A \code{\linkS4class{NeuroSpace}} object defining the 3D space dimensions,
spacing, and orientation.}

\item{indices}{An integer vector containing the linear indices of the voxels to include
in the lookup volume. These should be 1-based indices within the range of the space.}
}
\value{
An object of class \code{\linkS4class{IndexLookupVol}} containing:
\itemize{
  \item A mapping between linear indices and sparse positions
  \item The original space information
  \item The subset of included voxel indices
}
}
\description{
A three-dimensional brain image class that serves as a map between 1D grid indices
and a table of values. This class is primarily used in conjunction with the
\code{\linkS4class{SparseNeuroVec}} class to efficiently represent and access
sparse neuroimaging data.

The \code{IndexLookupVol} class provides efficient indexing and coordinate lookup
functionality for 3D neuroimaging data. It maintains a mapping between linear indices
and 3D coordinates, optimizing memory usage and access speed for sparse volumes.

Creates an \code{\linkS4class{IndexLookupVol}} object, which provides efficient
bidirectional mapping between linear indices and 3D coordinates in a neuroimaging
volume. This is particularly useful for working with masked or sparse brain volumes.
}
\details{
The IndexLookupVol class extends \code{\linkS4class{NeuroVol}} and provides a
mechanism for efficient lookup and mapping of sparse 3D neuroimaging data. It
stores only the indices of non-zero voxels and their corresponding mappings,
allowing for memory-efficient representation of large, sparse brain images.

Create an IndexLookupVol Object
}
\section{Slots}{

\describe{
\item{\code{space}}{A \code{\linkS4class{NeuroSpace}} object representing the 3D space of the brain image.}

\item{\code{indices}}{An integer vector containing the 1D indices of the non-zero voxels in the grid.}

\item{\code{map}}{An integer vector containing the mapping between the 1D indices and the table of values.}
}}

\section{Methods}{

This class inherits methods from \code{\linkS4class{NeuroVol}}. Additional
methods specific to index lookup and mapping operations may be available.
}

\section{Implementation Details}{

The class uses an integer mapping array for O(1) lookups between linear indices
and their corresponding positions in the sparse representation.
}

\examples{
# Create a NeuroSpace object
space <- NeuroSpace(dim = c(2L, 2L, 2L), origin = c(0, 0, 0), spacing = c(1, 1, 1))

# Create a 3D mask
mask <- array(c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE), dim = c(2, 2, 2))

# Create indices and map for the IndexLookupVol
indices <- which(mask)
map <- seq_along(indices)

# Create an IndexLookupVol object
ilv <- IndexLookupVol(space = space, indices = as.integer(indices))

# Access the indices
print(ilv@indices)

# Access the map
print(ilv@map)


# Create a 64x64x64 space
space <- NeuroSpace(c(64, 64, 64), c(1, 1, 1), c(0, 0, 0))

# Create a lookup volume with random indices
indices <- sample(1:262144, 10000)  # Select 10000 random voxels
ilv <- IndexLookupVol(space, indices)

# Look up coordinates for specific indices
coords <- coords(ilv, indices[1:10])


}
\seealso{
\code{\link{SparseNeuroVec-class}} for the primary class that utilizes IndexLookupVol.
\code{\link{NeuroVol-class}} for the base volumetric image class.

\code{\link{IndexLookupVol}} for creating instances of this class

\code{\link{coords}} for coordinate lookup,
\code{\link{lookup}} for index mapping,
\code{\linkS4class{NeuroSpace}} for space representation
}
