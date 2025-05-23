% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R, R/roi.R
\name{voxels}
\alias{voxels}
\alias{voxels,Kernel-method}
\title{extract voxel coordinates}
\usage{
voxels(x, ...)

\S4method{voxels}{Kernel}(x, center_voxel = NULL)
}
\arguments{
\item{x}{the object to extract voxels from}

\item{...}{additional arguments to function}

\item{center_voxel}{the absolute location of the center of the voxel, default is (0,0,0)}
}
\value{
A \code{matrix} or \code{vector} representing voxel coordinates from \code{x}.
}
\description{
extract voxel coordinates
}
\examples{
# Create a 3D kernel with dimensions 3x3x3 and voxel size 1x1x1
kern <- Kernel(kerndim = c(3,3,3), vdim = c(1,1,1))

# Get voxel coordinates centered at origin (0,0,0)
vox <- voxels(kern)
# Returns a matrix where each row is a voxel coordinate
# relative to the kernel center

# Get voxel coordinates centered at specific point (5,5,5)
vox_centered <- voxels(kern, center_voxel = c(5,5,5))
# Returns coordinates shifted to be centered at (5,5,5)

}
