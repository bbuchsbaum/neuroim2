# DenseNeuroVol Class

Represents a three-dimensional brain image backed by a dense array. This
class combines the spatial properties of
[`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
with the data storage capabilities of an array.

Construct a `DenseNeuroVol` instance

## Usage

``` r
DenseNeuroVol(data, space, label = "", indices = NULL)
```

## Arguments

- data:

  a three-dimensional `array`

- space:

  an instance of class
  [`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md)

- label:

  a `character` string

- indices:

  an optional 1-d index vector

## Value

`DenseNeuroVol` instance

## Details

DenseNeuroVol objects are used for 3D brain images where most or all
voxels contain meaningful data. They provide efficient access to
individual voxel values and are suitable for operations that require
frequent random access to voxel data.

## See also

[`NeuroVol-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md),
[`SparseNeuroVol-class`](https://bbuchsbaum.github.io/neuroim2/reference/SparseNeuroVol-class.md)

## Examples

``` r
# Create a simple 3D brain volume
vol_data <- array(rnorm(64*64*64), c(64, 64, 64))
vol_space <- NeuroSpace(dim=c(64L, 64L, 64L), origin=c(0, 0, 0), spacing=c(1, 1, 1))
brain_vol <- new("DenseNeuroVol", .Data=vol_data, space=vol_space)
```
