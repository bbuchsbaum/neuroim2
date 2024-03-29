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

\item{clusters}{a vector of clusters ids with length equal to number of nonzero voxels in mask \code{mask}}

\item{label_map}{an optional \code{list} that maps from cluster id to a cluster label, e.g. (1 -> "FFA", 2 -> "PPA")}

\item{label}{an optional \code{character} string used to label of the volume}
}
\value{
\code{\linkS4class{ClusteredNeuroVol}} instance
}
\description{
A three-dimensional brain image divided into N disjoint partitions or clusters.

Construct a \code{\linkS4class{ClusteredNeuroVol}} instance
}
\details{
The use case of \code{ClusteredNeuroVol} is to store volumetric data that has been clustered into discrete sets of voxels,
each of which has an associated id. For example, this class can be used to represent parcellated neuroimaging volumes.
}
\section{Slots}{

\describe{
\item{\code{mask}}{A \code{\linkS4class{LogicalNeuroVol}} object representing the logical mask indicating the spatial domain of the set of clusters.}

\item{\code{clusters}}{An integer vector representing the cluster number for each voxel in the mask.}

\item{\code{label_map}}{A named list where each element represents a cluster and its name.}

\item{\code{cluster_map}}{An \code{environment} object that maps from cluster id to the set of 1D spatial indices belonging to that cluster.}
}}

\section{Methods}{

  The \code{\linkS4class{ClusteredNeuroVol}} class inherits methods from the \code{\linkS4class{SparseNeuroVol}} class.
}

\section{Usage}{

  The \code{\linkS4class{ClusteredNeuroVol}} class is useful for representing a brain image with clustered regions, where each region is a disjoint partition.
}

\examples{

bspace <- NeuroSpace(c(16,16,16), spacing=c(1,1,1))
grid <- index_to_grid(bspace, 1:(16*16*16))
kres <- kmeans(grid, centers=10)
mask <- NeuroVol(rep(1, 16^3),bspace)
clusvol <- ClusteredNeuroVol(mask, kres$cluster)
}
\seealso{
\code{\linkS4class{SparseNeuroVol}}, \code{\linkS4class{LogicalNeuroVol}}
}
