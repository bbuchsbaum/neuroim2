% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_class.R, R/clustervol.R
\docType{class}
\name{ClusteredNeuroVol-class}
\alias{ClusteredNeuroVol-class}
\alias{ClusteredNeuroVol}
\title{ClusteredNeuroVol Class}
\usage{
ClusteredNeuroVol(mask, clusters, label_map = NULL, label = "")
}
\arguments{
\item{mask}{an instance of class \code{\linkS4class{LogicalNeuroVol}}}

\item{clusters}{a vector of clusters ids with length equal to number of nonzero
voxels in mask \code{mask}}

\item{label_map}{an optional \code{list} that maps from cluster id to a cluster
label, e.g. (1 -> "FFA", 2 -> "PPA")}

\item{label}{an optional \code{character} string used to label of the volume}
}
\value{
\code{\linkS4class{ClusteredNeuroVol}} instance
}
\description{
This class represents a three-dimensional brain image divided into N disjoint
partitions or clusters. It extends the \code{\linkS4class{SparseNeuroVol}} class
to provide efficient storage and manipulation of clustered neuroimaging data.

Construct a \code{\linkS4class{ClusteredNeuroVol}} instance
}
\details{
The ClusteredNeuroVol class is designed for efficient representation and
manipulation of brain images with distinct, non-overlapping regions or clusters.
It combines the memory efficiency of sparse representations with additional
structures for managing cluster information.

The use case of \code{ClusteredNeuroVol} is to store volumetric data that has been clustered into discrete sets of voxels,
each of which has an associated id. For example, this class can be used to represent parcellated neuroimaging volumes.
}
\section{Slots}{

\describe{
\item{\code{mask}}{A \code{\linkS4class{LogicalNeuroVol}} object representing the logical
mask indicating the spatial domain of the set of clusters.}

\item{\code{clusters}}{An integer vector representing the cluster number for each voxel
in the mask.}

\item{\code{label_map}}{A named list where each element represents a cluster and its name.}

\item{\code{cluster_map}}{An \code{environment} object that maps from cluster id to the
set of 1D spatial indices belonging to that cluster.}
}}

\section{Methods}{

This class inherits methods from the \code{\linkS4class{SparseNeuroVol}} class.
Additional methods specific to cluster operations may be available.
}

\section{Usage}{

ClusteredNeuroVol objects are particularly useful for:
\itemize{
  \item Representing parcellated brain images
  \item Storing results of clustering algorithms applied to neuroimaging data
  \item Efficient manipulation and analysis of region-based neuroimaging data
}
}

\examples{

# Create a simple clustered brain volume
dim <- c(10L, 10L, 10L)
mask_data <- array(rep(c(TRUE, FALSE), 500), dim)
mask <- new("LogicalNeuroVol", .Data = mask_data,
            space = NeuroSpace(dim = dim, origin = c(0,0,0), spacing = c(1,1,1)))

clusters <- as.integer(runif(sum(mask_data)) * 5)+1
label_map <- list("Cluster1" = 1, "Cluster2" = 2, "Cluster3" = 3,
                  "Cluster4" = 4, "Cluster5" = 5)

cluster_map <- list()
for (i in 1:5) {
  cluster_map[[as.character(i)]] <- which(clusters == i)
}

clustered_vol <- ClusteredNeuroVol(
                     mask = mask,
                     clusters = clusters,
                     label_map = label_map)



# Create a simple space and volume
space <- NeuroSpace(c(16, 16, 16), spacing = c(1, 1, 1))
vol_data <- array(rnorm(16^3), dim = c(16, 16, 16))
vol <- NeuroVol(vol_data, space)

# Create a binary mask (e.g., values > 0)
mask_data <- vol_data > 0
mask_vol <- LogicalNeuroVol(mask_data, space)

# Get coordinates of masked voxels
mask_idx <- which(mask_data)
coords <- index_to_coord(mask_vol, mask_idx)

# Cluster the coordinates into 10 groups
set.seed(123)  # for reproducibility
kmeans_result <- kmeans(coords, centers = 10)

# Create the clustered volume
clustered_vol <- ClusteredNeuroVol(mask_vol, kmeans_result$cluster)

# Print information about the clusters
print(clustered_vol)
}
\seealso{
\code{\link{SparseNeuroVol-class}} for the parent sparse volume class.
\code{\link{LogicalNeuroVol-class}} for the mask representation.
}
