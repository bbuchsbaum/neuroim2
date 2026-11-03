# SparseNeuroVol Class

This class represents a three-dimensional brain image using a sparse
data representation. It is particularly useful for large brain images
with a high proportion of zero or missing values, offering efficient
storage and processing.

Construct a `SparseNeuroVol` instance

## Usage

``` r
SparseNeuroVol(data, space, indices = NULL, label = "")
```

## Arguments

- data:

  a numeric vector or ROIVol

- space:

  an instance of class
  [`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md)

- indices:

  a index vector indicating the 1-d coordinates of the data values

- label:

  a `character` string

## Value

`SparseNeuroVol` instance

## Details

The SparseNeuroVol class extends the
[`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
class and implements the ArrayLike3D interface. It uses a `sparseVector`
from the Matrix package to store the image data, which allows for
memory-efficient representation of sparse 3D neuroimaging data.

Image data is backed by
[`Matrix::sparseVector`](https://rdrr.io/pkg/Matrix/man/sparseVector.html).

## Slots

- `data`:

  A `sparseVector` object from the Matrix package, storing the image
  volume data in a sparse format.

## References

Bates, D., & Maechler, M. (2019). Matrix: Sparse and Dense Matrix
Classes and Methods. R package version 1.2-18.
https://CRAN.R-project.org/package=Matrix

## See also

[`NeuroVol-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
for the base volumetric image class.
[`DenseNeuroVol-class`](https://bbuchsbaum.github.io/neuroim2/reference/DenseNeuroVol-class.md)
for a dense representation of 3D brain images.

## Examples

``` r

# Create a sparse 3D brain image
dim <- c(64L, 64L, 64L)
space <- NeuroSpace(dim = dim, origin = c(0, 0, 0), spacing = c(1, 1, 1))
sparse_data <- Matrix::sparseVector(x = c(1, 2, 3),
                                    i = c(100, 1000, 10000),
                                    length = prod(dim))
sparse_vol <- new("SparseNeuroVol", space = space, data = sparse_data)
sparse_vol[1000] == 1
#> [1] FALSE

data <- 1:10
indices <- seq(1,1000, length.out=10)
bspace <- NeuroSpace(c(64,64,64), spacing=c(1,1,1))
sparsevol <- SparseNeuroVol(data,bspace,indices=indices)
densevol <- NeuroVol(data,bspace,indices=indices)
sum(sparsevol) == sum(densevol)
#> [1] TRUE
```
