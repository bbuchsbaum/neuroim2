# Logical Negation for Neuroimaging Volumes

Logical negation (`!`) for neuroimaging volumes. Returns a
[`LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.md)
where non-zero voxels become `FALSE` and zero voxels become `TRUE`.

## Usage

``` r
# S4 method for class 'DenseNeuroVol'
!x

# S4 method for class 'SparseNeuroVol'
!x
```

## Arguments

- x:

  A neuroimaging volume object.

## Value

A
[`LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.md).

## Examples

``` r
sp <- NeuroSpace(c(5L, 5L, 5L), c(1, 1, 1))
mask <- LogicalNeuroVol(array(sample(c(TRUE, FALSE), 125, replace = TRUE),
                              c(5, 5, 5)), sp)
inv <- !mask
```
