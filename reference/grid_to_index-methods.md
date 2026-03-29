# Generic function to convert N-dimensional grid coordinates to 1D indices

Converts 2D grid coordinates to linear indices for a `NeuroSlice`
object.

## Usage

``` r
grid_to_index(x, coords)

# S4 method for class 'NeuroSlice,matrix'
grid_to_index(x, coords)

# S4 method for class 'NeuroSlice,numeric'
grid_to_index(x, coords)

# S4 method for class 'NeuroSpace,matrix'
grid_to_index(x, coords)

# S4 method for class 'NeuroSpace,numeric'
grid_to_index(x, coords)

# S4 method for class 'NeuroVol,matrix'
grid_to_index(x, coords)

# S4 method for class 'NeuroVol,numeric'
grid_to_index(x, coords)
```

## Arguments

- x:

  A `NeuroSlice` object

- coords:

  Either a numeric vector of length 2 or a matrix with 2 columns,
  representing (x,y) coordinates in the slice grid

## Value

An integer `vector` of 1D indices corresponding to `coords`.

## Details

Convert Grid Coordinates to Linear Indices

## See also

[`index_to_grid`](https://bbuchsbaum.github.io/neuroim2/reference/index_to_grid-methods.md)
for the inverse operation

## Examples

``` r
# Create a 2D space (10x10)
space_2d <- NeuroSpace(c(10,10), c(1,1))

# Convert 2D grid coordinates to linear indices
coords_2d <- matrix(c(1,1, 2,2), ncol=2, byrow=TRUE)
idx_2d <- grid_to_index(space_2d, coords_2d)
# First coordinate (1,1) maps to index 1
# Second coordinate (2,2) maps to index 12 (= 2 + (2-1)*10)

# Create a 3D space (10x10x10)
space_3d <- NeuroSpace(c(10,10,10), c(1,1,1))

# Convert 3D grid coordinates to linear indices
coords_3d <- matrix(c(1,1,1, 2,2,2), ncol=3, byrow=TRUE)
idx_3d <- grid_to_index(space_3d, coords_3d)

# Single coordinate can also be converted
idx <- grid_to_index(space_3d, c(1,1,1))

slice_space <- NeuroSpace(c(10, 10))
slice_data <- matrix(1:100, 10, 10)
slice <- NeuroSlice(slice_data, slice_space)

# Convert single coordinate
idx <- grid_to_index(slice, c(5, 5))

# Convert multiple coordinates
coords <- matrix(c(1,1, 2,2, 3,3), ncol=2, byrow=TRUE)
indices <- grid_to_index(slice, coords)
```
