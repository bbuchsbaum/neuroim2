# Number of Clusters

This function returns the number of clusters in a ClusteredNeuroVol
object.

## Usage

``` r
num_clusters(x)

# S4 method for class 'ClusteredNeuroVec'
num_clusters(x)

# S4 method for class 'ClusteredNeuroVol'
num_clusters(x)
```

## Arguments

- x:

  A ClusteredNeuroVol object.

## Value

An `integer` representing the number of clusters in `x`.

An integer representing the number of clusters in the input object.

## Examples

``` r
# Create a simple 3D volume and mask
space <- NeuroSpace(c(16, 16, 16), spacing = c(1, 1, 1))
vol_data <- array(rnorm(16^3), dim = c(16, 16, 16))
mask_vol <- LogicalNeuroVol(vol_data > 0, space)

# Get coordinates of masked voxels for clustering
mask_idx <- which(mask_vol)
coords <- index_to_coord(mask_vol, mask_idx)

# Cluster the coordinates into 10 groups using k-means
set.seed(123)  # for reproducibility
kmeans_result <- kmeans(coords, centers = 10)

# Create a clustered volume
clustered_vol <- ClusteredNeuroVol(mask_vol, kmeans_result$cluster)

# Get the number of clusters
n_clusters <- num_clusters(clustered_vol)
n_clusters == 10
#> [1] TRUE
```
