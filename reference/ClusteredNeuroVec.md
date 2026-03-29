# ClusteredNeuroVec: Cluster-aware 4D neuroimaging data

\`ClusteredNeuroVec\` creates a 4D array-like object where voxels are
grouped into clusters, with one time-series per cluster. All voxels
within a cluster share the same time-series, making it ideal for
parcellated analyses (e.g., Schaefer-Yeo parcellations).

## Usage

``` r
ClusteredNeuroVec(x, cvol, FUN = mean, weights = NULL, label = "")
```

## Arguments

- x:

  Either a \`NeuroVec\` object to be reduced by clusters, or a
  pre-computed numeric matrix of cluster time-series (T x K, where
  T=time points, K=clusters)

- cvol:

  A \`ClusteredNeuroVol\` object defining the cluster assignments

- FUN:

  Reduction function to aggregate voxels within clusters (default:
  mean). Common choices include `mean`, `median`, or custom functions.

- weights:

  Optional numeric vector of per-voxel weights for weighted aggregation.
  Must have length equal to the number of non-zero voxels in the mask.

- label:

  Optional character label for the object (default: "")

## Value

A `ClusteredNeuroVec` object containing:

- cvol:

  The input ClusteredNeuroVol defining cluster structure

- ts:

  A TxK matrix of cluster time-series (T=timepoints, K=clusters)

- cl_map:

  Integer vector mapping linear voxel indices to cluster IDs

- label:

  Character label for the object

## Details

This class implements array-like 4D access while storing data
efficiently as a TxK matrix instead of the full voxel x time
representation. Each cluster's time-series is computed by applying the
aggregation function (`FUN`) to all voxels within that cluster.

The object supports standard NeuroVec operations:

- Indexing: `x[,,,t]` to extract 3D volumes at time t

- Series extraction: `series(x, i, j, k)` for time-series at voxel
  (i,j,k)

- Matrix conversion: `as.matrix(x)` to get the TxK cluster matrix

Single-voxel clusters are handled efficiently without aggregation
overhead.

## See also

[`ClusteredNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/ClusteredNeuroVol-class.md)
for creating cluster assignments,
[`cluster_searchlight_series`](https://bbuchsbaum.github.io/neuroim2/reference/cluster_searchlight_series.md)
for cluster-based searchlight analysis,
[`series`](https://bbuchsbaum.github.io/neuroim2/reference/series-methods.md)
for extracting time-series

## Examples

``` r
# Create synthetic 4D data (10x10x10 volume, 20 timepoints)
sp4 <- NeuroSpace(c(10,10,10,20), c(1,1,1))
arr <- array(rnorm(10*10*10*20), dim=c(10,10,10,20))
vec <- NeuroVec(arr, sp4)

# Create a mask covering the central region
sp3 <- NeuroSpace(c(10,10,10), c(1,1,1))
mask_arr <- array(FALSE, dim=c(10,10,10))
mask_arr[3:8, 3:8, 3:8] <- TRUE
mask <- LogicalNeuroVol(mask_arr, sp3)

# Assign voxels to 5 random clusters
n_voxels <- sum(mask_arr)
clusters <- sample(1:5, n_voxels, replace=TRUE)
cvol <- ClusteredNeuroVol(mask, clusters)

# Create clustered representation
cv <- ClusteredNeuroVec(vec, cvol)

# Access like a regular NeuroVec
vol_t1 <- cv[,,,1]  # 3D volume at time 1
ts <- series(cv, 5, 5, 5)  # time-series at voxel (5,5,5)

# Get cluster time-series matrix
cluster_matrix <- as.matrix(cv)  # T x K matrix
dim(cluster_matrix)  # 20 x 5
#> [1] 20  5
```
