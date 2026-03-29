# Cut an object into a list of spatial or spatiotemporal clusters

This function cuts an object into a list of sub-objects based on a
vector of cluster indices. The resulting list contains each of the
clusters as separate objects.

These methods split a NeuroVec object into multiple ROIVec objects based
on cluster assignments.

## Usage

``` r
split_clusters(x, clusters, ...)

# S4 method for class 'NeuroVec,ClusteredNeuroVol'
split_clusters(x, clusters, ...)

# S4 method for class 'NeuroVec,integer'
split_clusters(x, clusters, ...)

# S4 method for class 'NeuroVol,ClusteredNeuroVol'
split_clusters(x, clusters)

# S4 method for class 'NeuroVol,integer'
split_clusters(x, clusters)

# S4 method for class 'NeuroVol,numeric'
split_clusters(x, clusters)

# S4 method for class 'ClusteredNeuroVol,missing'
split_clusters(x, clusters)

# S4 method for class 'NeuroVec,integer'
split_clusters(x, clusters, ...)

# S4 method for class 'NeuroVec,numeric'
split_clusters(x, clusters, ...)

# S4 method for class 'NeuroVec,ClusteredNeuroVol'
split_clusters(x, clusters, ...)
```

## Arguments

- x:

  A NeuroVec object to be split.

- clusters:

  Either a ClusteredNeuroVol object or an integer vector of cluster
  assignments.

- ...:

  Additional arguments to be passed to methods.

## Value

A `list` of sub-objects, where each sub-object corresponds to a unique
cluster index.

A deflist (lazy-loading list) of ROIVec objects, where each element
corresponds to a cluster.

## Details

There are two methods for splitting clusters:

- Using a ClusteredNeuroVol object: This method uses the pre-defined
  clusters in the ClusteredNeuroVol object.

- Using an integer vector: This method allows for custom cluster
  assignments.

methods return a deflist, which is a lazy-loading list of ROIVec
objects.

## See also

[`NeuroVec-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md),
[`ClusteredNeuroVol-class`](https://bbuchsbaum.github.io/neuroim2/reference/ClusteredNeuroVol-class.md),
[`ROIVec-class`](https://bbuchsbaum.github.io/neuroim2/reference/ROIVec-class.md)

## Examples

``` r
# \donttest{
  # Create a synthetic 3D volume and its NeuroSpace
  space <- NeuroSpace(c(10, 10, 10,4))
  vol_data <- array(rnorm(10 * 10 * 10 * 4), dim = c(10, 10, 10,4))
  neuro_vec <- NeuroVec(vol_data, space)

  # Create a binary mask (e.g., select voxels with values > 0)
  mask_data <- as.logical(neuro_vec[[1]] > .5)
  mask_vol <- LogicalNeuroVol(mask_data, NeuroSpace(c(10, 10, 10)))

  # Extract indices and coordinates for the masked voxels
  mask_idx <- which(mask_data)
  coords <- index_to_coord(mask_vol, mask_idx)

  # Perform k-means clustering on the coordinates (e.g., 3 clusters)
  set.seed(123)  # for reproducibility
  k_res <- kmeans(coords, centers = 3)

  # Create a ClusteredNeuroVol using the mask and k-means cluster assignments
  clustered_vol <- ClusteredNeuroVol(mask_vol, k_res$cluster)

  # Split the NeuroVec by clusters using the ClusteredNeuroVol method
  split_result_clust <- split_clusters(neuro_vec, clustered_vol)

  # Calculate and print the mean value for each cluster
  means_clust <- sapply(split_result_clust, function(x) mean(values(x)))
  print(means_clust)
#> [1] 0.3150001 0.2943747 0.3140765

  # Alternatively, create an integer vector of cluster assignments:
  cluster_assignments <- numeric(prod(dim(space)[1:3]))
  cluster_assignments[mask_idx] <- k_res$cluster
  split_result_int <- split_clusters(neuro_vec, as.integer(cluster_assignments))

  # Verify that both splitting methods yield the same cluster means
  means_int <- sapply(split_result_int, function(x) mean(values(x)))
  print(all.equal(sort(means_clust), sort(means_int)))
#> [1] TRUE
# }


# Create a simple example space and data
space <- NeuroSpace(c(10, 10, 10,4))
data <- array(rnorm(1000*4), dim = c(10, 10, 10,4))
vec <- NeuroVec(data, space)

# Create a mask for clustering (e.g., values > 0)
mask <- vec[,,,1] > 0
mask_vol <- LogicalNeuroVol(as.array(mask), NeuroSpace(c(10, 10, 10)))

# Get coordinates of masked voxels for clustering
mask_idx <- which(mask)
coords <- index_to_coord(mask_vol, mask_idx)

# Perform clustering on the coordinates (3 clusters for example)
set.seed(123) # for reproducibility
kmeans_result <- kmeans(coords, centers = 3)

# Create a ClusteredNeuroVol
clustered_vol <- ClusteredNeuroVol(mask_vol, kmeans_result$cluster)

# Split the NeuroVec by clusters
split_result <- split_clusters(vec, clustered_vol)

# Calculate mean value for each cluster
cluster_means <- sapply(split_result, function(x) mean(values(x)))
print(cluster_means)
#> [1] 0.1735834 0.1777742 0.2143060

# Alternative: using integer cluster assignments
cluster_indices <- numeric(prod(dim(space)[1:3]))
cluster_indices[mask_idx] <- kmeans_result$cluster
split_result2 <- split_clusters(vec, as.integer(cluster_indices))

# Verify both methods give same results
cluster_means2 <- sapply(split_result2, function(x) mean(values(x)))
print(all.equal(sort(cluster_means), sort(cluster_means2)))
#> [1] TRUE
```
