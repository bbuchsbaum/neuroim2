# ClusteredNeuroVec: Parcel-based 4D Analysis

## Overview

`ClusteredNeuroVec` provides an efficient representation for parcellated
4D neuroimaging data where voxels are grouped into clusters or parcels.
Instead of storing a time-series for every voxel, it stores one
time-series per cluster, making it ideal for analyses using brain
atlases like Schaefer-Yeo parcellations.

## Why ClusteredNeuroVec?

Traditional neuroimaging analyses often involve: - Reducing voxel-level
data to parcel/ROI averages - Working with brain atlases that group
voxels into regions - Performing searchlight analyses at the parcel
level rather than voxel level

`ClusteredNeuroVec` makes these workflows more efficient while
maintaining compatibility with standard `NeuroVec` operations.

## Creating a ClusteredNeuroVec

### From scratch with synthetic data

``` r
# Create a simple 3D space with mask
space <- NeuroSpace(c(10, 10, 10), spacing = c(2, 2, 2))
mask_data <- array(TRUE, c(10, 10, 10))
mask_data[1:3, 1:3, 1:3] <- FALSE  # exclude corner
mask <- LogicalNeuroVol(mask_data, space)

# Create cluster assignments (e.g., 5 random clusters)
n_masked <- sum(mask_data)
cluster_ids <- sample(1:5, n_masked, replace = TRUE)
cvol <- ClusteredNeuroVol(mask, cluster_ids)

# Create synthetic 4D data
vec_space <- NeuroSpace(c(10, 10, 10, 20), spacing = c(2, 2, 2))
vec_data <- array(rnorm(10 * 10 * 10 * 20), dim = c(10, 10, 10, 20))
vec <- NeuroVec(vec_data, vec_space)

# Create ClusteredNeuroVec
cv <- ClusteredNeuroVec(vec, cvol)
print(cv)
#> <ClusteredNeuroVec> [5 clusters x 20 timepoints] 
#> ── Spatial ───────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 10 x 10 x 10
#>   Spacing       : 2 x 2 x 2 mm
#> ── Cluster Data ──────────────────────────────────────────────────────────────── 
#>   Clusters      : 5
#>   Time Points   : 20
#>   Cluster Sizes : min=180, median=199, max=202
#>   Size          : 42.9 Kb
```

### Key properties

``` r
# Dimensions: still 4D (x, y, z, time)
dim(cv)
#> [1] 10 10 10 20

# Number of clusters
num_clusters(cv)
#> [1] 5

# Access cluster time-series matrix (T x K)
ts_matrix <- as.matrix(cv, by = "cluster")
dim(ts_matrix)  # 20 time points x 5 clusters
#> [1] 20  5
```

## Array-like access

`ClusteredNeuroVec` behaves like a regular 4D array:

``` r
# Extract 3D volume at time point 1
vol_t1 <- cv[,,,1]
dim(vol_t1)  # 10 x 10 x 10
#> [1] 10 10 10

# All voxels in the same cluster have the same value
# (they share the cluster's mean time-series)

# Get time-series at a specific voxel
ts <- series(cv, 5, 5, 5)
length(ts)  # 20 time points
#> [1] 20
```

## Cluster searchlight analysis

Perform searchlight analysis at the cluster level using centroid
distances:

``` r
# K-nearest neighbor searchlight (3 nearest clusters)
windows_knn <- cluster_searchlight_series(cv, k = 3)
length(windows_knn)  # One window per cluster
#> [1] 5

# Look at first window
win1 <- windows_knn[[1]]
dim(values(win1))  # 20 time points x 3 neighbors
#> [1] 20  3

# Radius-based searchlight (e.g., 15mm radius)
windows_radius <- cluster_searchlight_series(cv, radius = 15)
```

## Real-world example: Schaefer parcellation

``` r
# Load fMRI data
fmri_data <- read_vec("subject01_task.nii.gz")

# Load Schaefer atlas (example with 400 parcels)
atlas <- read_vol("Schaefer2018_400Parcels_7Networks.nii.gz")
mask <- atlas > 0

# Create ClusteredNeuroVol from atlas
cvol <- ClusteredNeuroVol(mask, as.integer(atlas[mask]))

# Create parcellated representation
cv <- ClusteredNeuroVec(fmri_data, cvol)

# Now you have 400 time-series (one per parcel) instead of ~200,000 voxels
parcels <- as.matrix(cv, by = "cluster")
dim(parcels)  # T x 400

# Perform connectivity analysis at parcel level
cor_matrix <- cor(parcels)
dim(cor_matrix)  # 400 x 400
```

## Integration with existing workflows

`ClusteredNeuroVec` integrates seamlessly with existing neuroim2
functions:

``` r
# Use with split_reduce for custom aggregation
# (ClusteredNeuroVec already uses this internally)

# Scale time-series within each cluster
# (if scale_series is implemented for ClusteredNeuroVec)
# cv_scaled <- scale_series(cv, center = TRUE, scale = TRUE)

# Get cluster centroids for visualization
centers <- centroids(cv)
head(centers)  # x, y, z coordinates
#>          [,1]     [,2]     [,3]
#> [1,] 5.783920 5.356784 5.522613
#> [2,] 5.495050 5.544554 5.579208
#> [3,] 5.406250 5.588542 5.604167
#> [4,] 5.683333 5.600000 5.638889
#> [5,] 5.620000 5.895000 5.645000
```

## Performance benefits

By storing only K time-series instead of N voxels: - Memory usage: O(K ×
T) instead of O(N × T) - Searchlight operations: O(K²) instead of
O(N²) - Typical reduction: 100-1000x fewer time-series

For a typical fMRI dataset: - Voxel-level: ~200,000 voxels × 500
timepoints = 100M values - Parcel-level: 400 parcels × 500 timepoints =
200K values

## Summary

`ClusteredNeuroVec` provides: - Efficient storage for parcellated 4D
data - Full array-like access semantics - Cluster-aware searchlight
operations - Seamless integration with existing neuroim2 workflows

It’s ideal for: - Atlas-based analyses - Connectivity studies -
Parcellated machine learning - Any workflow that aggregates voxels to
regions
