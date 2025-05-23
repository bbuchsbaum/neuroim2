% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R, R/neurospace.R
\name{grid_to_grid}
\alias{grid_to_grid}
\alias{grid_to_grid,NeuroSpace,matrix-method}
\alias{grid_to_grid,matrix,matrix-method}
\title{Generic function to convert voxel coordinates in the reference space (LPI) to native array space.}
\usage{
grid_to_grid(x, vox)

\S4method{grid_to_grid}{NeuroSpace,matrix}(x, vox)

\S4method{grid_to_grid}{matrix,matrix}(x, vox)
}
\arguments{
\item{x}{the object}

\item{vox}{a matrix of LPI voxel coordinates}
}
\value{
A numeric \code{matrix} of native voxel coordinates.
}
\description{
Generic function to convert voxel coordinates in the reference space (LPI) to native array space.
}
\examples{
# Create a simple 3D volume in LPI orientation
space <- NeuroSpace(c(10,10,10), c(2,2,2))

# Create a reoriented space in RAS orientation
space_ras <- reorient(space, c("R", "A", "S"))

# Convert coordinates between orientations
voxel_coords <- t(matrix(c(1,1,1)))
new_coords <- grid_to_grid(space_ras, voxel_coords)
print(new_coords)
}
