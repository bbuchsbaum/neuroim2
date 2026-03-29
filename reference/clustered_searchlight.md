# Create a clustered searchlight iterator

This function generates a searchlight iterator that iterates over
successive spatial clusters in an image volume. It allows for the
exploration of spatially clustered regions within the provided mask by
using either a pre-defined clustered volume or performing k-means
clustering to generate the clusters.

## Usage

``` r
clustered_searchlight(mask, cvol = NULL, csize = NULL)
```

## Arguments

- mask:

  A
  [`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
  object representing the brain mask.

- cvol:

  An optional `ClusteredNeuroVol` instance representing pre-defined
  clusters within the mask. If provided, the 'csize' parameter is
  ignored.

- csize:

  An optional integer specifying the number of clusters to be generated
  using k-means clustering (ignored if `cvol` is provided).

## Value

A `deferred_list` object containing `ROIVol` objects, each representing
a clustered region within the image volume.

## Examples

``` r
# Load an example brain mask
mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))

# Generate a clustered searchlight iterator with 5 clusters
clust_searchlight <- clustered_searchlight(mask, csize = 5)

```
