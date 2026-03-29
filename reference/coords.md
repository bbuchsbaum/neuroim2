# Extract coordinates from an object

This function extracts the coordinates from an input object.

## Usage

``` r
coords(x, ...)
```

## Arguments

- x:

  The object to extract coordinates from.

- ...:

  Additional arguments (not used in the generic function).

## Value

A numeric `matrix` or `vector` containing the coordinates of `x`.

## Examples

``` r
# Create a NeuroSpace object with 3mm voxels
space <- NeuroSpace(c(10,10,10), spacing=c(3,3,3))

# Create ROI coordinates in voxel space
coords <- matrix(c(1,1,1, 2,2,2), ncol=3, byrow=TRUE)
roi_coords <- ROICoords(coords)

# Get coordinates in voxel space
vox_coords <- coords(roi_coords)
# First coordinate is (1,1,1)

# Get coordinates
cds <- coords(roi_coords)
nrow(cds) == 2
#> [1] TRUE
```
