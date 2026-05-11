# Convert 1d indices to n-dimensional grid coordinates

Converts linear indices to 2D grid coordinates for a `NeuroSlice`
object.

## Usage

``` r
index_to_grid(x, idx)

# S4 method for class 'NeuroSlice,numeric'
index_to_grid(x, idx)

# S4 method for class 'NeuroSpace,numeric'
index_to_grid(x, idx)

# S4 method for class 'NeuroVec,index'
index_to_grid(x, idx)

# S4 method for class 'NeuroVec,integer'
index_to_grid(x, idx)

# S4 method for class 'NeuroVol,index'
index_to_grid(x, idx)

# S4 method for class 'NeuroVol,integer'
index_to_grid(x, idx)
```

## Arguments

- x:

  A `NeuroSlice` object

- idx:

  Integer vector of linear indices to convert

## Value

A numeric `matrix` of grid coordinates.

## Details

Convert Linear Indices to Grid Coordinates

## See also

[`grid_to_index`](https://bbuchsbaum.github.io/neuroim2/reference/grid_to_index-methods.md)
for the inverse operation

## Examples

``` r

 bvol <- NeuroVol(array(0, c(10,10,10)), NeuroSpace(c(10,10,10), c(1,1,1)))
 idx <- 1:10
 g <- index_to_grid(bvol, idx)
 bvol[g]
#>  [1] 0 0 0 0 0 0 0 0 0 0

slice_space <- NeuroSpace(c(10, 10))
slice_data <- matrix(1:100, 10, 10)
slice <- NeuroSlice(slice_data, slice_space)

# Convert single index
coords <- index_to_grid(slice, 55)

# Convert multiple indices
indices <- c(1, 25, 50, 75, 100)
coords_mat <- index_to_grid(slice, indices)
```
