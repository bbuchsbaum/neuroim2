# Create an instance of class [`ROIVec`](https://bbuchsbaum.github.io/neuroim2/reference/ROIVec-class.md)

This function constructs an instance of the ROIVec class, which
represents a region of interest (ROI) in a 4D volume. The class stores
the NeuroSpace object, voxel coordinates, and data values for the ROI.

## Usage

``` r
ROIVec(vspace, coords, data = matrix(1, nrow = 1, ncol = nrow(coords)))
```

## Arguments

- vspace:

  An instance of class `NeuroSpace` with four dimensions, which
  represents the dimensions, voxel spacing, and time points of the 4D
  volume.

- coords:

  A 3-column matrix of voxel coordinates for the region of interest.

- data:

  The matrix of data values associated with the region of interest, with
  each row representing a voxel and each column representing a time
  point. By default, it is a matrix with a number of rows equal to the
  number of rows in the \`coords\` matrix and a single column filled
  with ones.

## Value

An instance of class `ROIVec`, containing the NeuroSpace object, voxel
coordinates, and data values for the region of interest.

## Examples

``` r
# Create a NeuroSpace object
vspace <- NeuroSpace(dim = c(5, 5, 5, 10), spacing = c(1, 1, 1))

# Define voxel coordinates for the ROI
coords <- matrix(c(1, 2, 3, 2, 2, 2, 3, 3, 3), ncol = 3)

# Create a data matrix for the ROI
data <- matrix(rnorm(30), nrow = 10, ncol = 3)

# Create a ROIVec object
roi_vec <- ROIVec(vspace, coords, data)
```
