# NeuroVol: 3D Neuroimaging Volume Class

The `NeuroVol` class encapsulates 3D volumetric neuroimaging data. It
provides methods for accessing slices, performing spatial
transformations, and integrating with the spatial reference provided by
[`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md).

## Usage

``` r
NeuroVol(data, space, label = "", indices = NULL)
```

## Arguments

- data:

  A 3D array containing the volumetric data.

- space:

  An object of class
  [`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md)
  defining the spatial properties.

- label:

  A character string providing a label for the volume (default: "").

- indices:

  An optional vector of indices for sparse representation (default:
  NULL).

## Value

A `NeuroVol` object.

## Examples

``` r
bspace <- NeuroSpace(c(64,64,64), spacing=c(1,1,1))
dat <- array(rnorm(64*64*64), c(64,64,64))
bvol <- NeuroVol(dat,bspace, label="test")
```
