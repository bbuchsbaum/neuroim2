# Cluster-centroid searchlight over cluster time-series

Iterate over clusters by their centroids and, for each seed cluster,
return the time-series of the K nearest clusters (or those within a
radius). This enables searchlight analysis at the cluster level rather
than individual voxels.

## Usage

``` r
cluster_searchlight_series(x, cvol = NULL, k = 10L, radius = NULL, label = "")
```

## Arguments

- x:

  A \`ClusteredNeuroVec\` object or a \`NeuroVec\` plus \`cvol\`

- cvol:

  A \`ClusteredNeuroVol\` (required if \`x\` is a \`NeuroVec\`)

- k:

  Integer, number of nearest clusters including the seed (default: 10).
  Will be capped at the total number of clusters if specified value
  exceeds it.

- radius:

  Numeric distance in mm. If given, use all clusters within this radius
  instead of k-nearest neighbors. Cannot be used together with `k`.

- label:

  Optional character label for the returned windows

## Value

A list of `ROIVec` objects, one per cluster, where each ROIVec contains:

- values:

  A TxN matrix where T is the number of timepoints and N is the number
  of neighboring clusters (including the seed itself)

- coords:

  The centroid coordinates of the neighboring clusters

The seed cluster's time-series is always the first column in each
ROIVec.

## Details

The function creates a searchlight around each cluster's centroid,
selecting either:

- The k nearest clusters (when `k` is specified)

- All clusters within a given radius (when `radius` is specified)

This is particularly useful for cluster-level connectivity analyses or
when working with parcellated data where voxel-level searchlights would
be redundant.

## See also

[`ClusteredNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/ClusteredNeuroVec.md)
for creating clustered neuroimaging vectors,
[`searchlight`](https://bbuchsbaum.github.io/neuroim2/reference/searchlight.md)
for voxel-level searchlight analysis,
[`ROIVec`](https://bbuchsbaum.github.io/neuroim2/reference/ROIVec.md)
for the structure of returned windows

## Examples

``` r
# Create synthetic 4D data (8x8x8 volume, 10 timepoints)
sp4 <- NeuroSpace(c(8,8,8,10), c(1,1,1))
arr <- array(rnorm(8*8*8*10), dim=c(8,8,8,10))
vec <- NeuroVec(arr, sp4)

# Create a mask covering most of the volume
sp3 <- NeuroSpace(c(8,8,8), c(1,1,1))
mask_arr <- array(FALSE, dim=c(8,8,8))
mask_arr[2:7, 2:7, 2:7] <- TRUE
mask <- LogicalNeuroVol(mask_arr, sp3)

# Assign voxels to 10 clusters
n_voxels <- sum(mask_arr)
clusters <- sample(1:10, n_voxels, replace=TRUE)
cvol <- ClusteredNeuroVol(mask, clusters)

# Create clustered representation
cv <- ClusteredNeuroVec(vec, cvol)

# Get cluster searchlight with 3 nearest neighbors
windows <- cluster_searchlight_series(cv, k = 3)
length(windows)  # 10 windows (one per cluster)
#> [1] 10

# Check first window
roi1 <- windows[[1]]
dim(values(roi1))  # 10 x 3 (timepoints x neighbors)
#> [1] 10  3

# Use radius-based neighborhoods (5mm radius)
windows_radius <- cluster_searchlight_series(cv, radius = 5)
# Each window may have different number of neighbors
```
