# NeuroVecSource

This function constructs a NeuroVecSource object, which represents the
source of a four-dimensional brain image.

## Usage

``` r
NeuroVecSource(file_name, indices = NULL, mask = NULL)
```

## Arguments

- file_name:

  The name of the 4-dimensional image file.

- indices:

  An optional integer vector specifying the subset of volume indices to
  load. If not provided, all volumes will be loaded.

- mask:

  An optional logical array or
  [`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
  object defining the subset of voxels to load. If provided, a
  SparseNeuroVecSource object will be created.

## Value

An instance of the
[`NeuroVecSource`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVecSource-class.md)
class.

## Details

If a `mask` is supplied, it should be a
[`LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.md)
or
[`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
instance. If the latter, then the mask will be defined by nonzero
elements of the volume.
