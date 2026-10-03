# Convert ClusteredNeuroVol to DenseNeuroVol

This method converts a ClusteredNeuroVol into an equivalent
DenseNeuroVol object.

## Usage

``` r
# S4 method for class 'ClusteredNeuroVol,DenseNeuroVol'
coerce(from, to = "DenseNeuroVol", strict = TRUE)
```

## Arguments

- from:

  A
  [`ClusteredNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/ClusteredNeuroVol-class.md)
  object to be converted

- to:

  The target class, `"DenseNeuroVol"`.

- strict:

  Logical; whether to enforce the target class exactly.

## Value

A
[`DenseNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/DenseNeuroVol-class.md)
object

## Details

Convert a ClusteredNeuroVol Object to a DenseNeuroVol Object

## See also

[`ClusteredNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/ClusteredNeuroVol-class.md),
[`DenseNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/DenseNeuroVol-class.md)

## Examples

``` r

# Create a clustered volume
mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))
clusters <- rep(1:5, length.out=sum(mask))
cvol <- ClusteredNeuroVol(mask, clusters)

# Convert to DenseNeuroVol
dvol <- as(cvol, "DenseNeuroVol")
```
