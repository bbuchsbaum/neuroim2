# Partition an image into a set of disjoint clusters

This function partitions an image into a set of disjoint clusters using
k-means clustering.

## Usage

``` r
partition(x, k, ...)

# S4 method for class 'LogicalNeuroVol,integer'
partition(x, k)

# S4 method for class 'LogicalNeuroVol,numeric'
partition(x, k)

# S4 method for class 'DenseNeuroVol,numeric'
partition(x, k)
```

## Arguments

- x:

  the image to partition, represented as a 3D array.

- k:

  the number of clusters to form.

- ...:

  additional arguments passed to the kmeans function.

## Value

a 3D `array` where each voxel is assigned to a cluster.

## See also

[`kmeans`](https://rdrr.io/r/stats/kmeans.html)

## Examples

``` r
# Load an example 3D image
library(neuroim2)
img <- read_vol(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))

# Partition the image into 5 clusters using default options
clusters <- partition(img, 5)

```
