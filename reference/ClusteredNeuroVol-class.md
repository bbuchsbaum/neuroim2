# ClusteredNeuroVol Class

This class represents a three-dimensional brain image divided into N
disjoint partitions or clusters. It extends the
[`SparseNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/SparseNeuroVol-class.md)
class to provide efficient storage and manipulation of clustered
neuroimaging data.

Construct a `ClusteredNeuroVol` instance

## Usage

``` r
ClusteredNeuroVol(mask, clusters, label_map = NULL, label = "")
```

## Arguments

- mask:

  an instance of class
  [`LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.md)

- clusters:

  a vector of clusters ids with length equal to number of nonzero voxels
  in mask `mask`

- label_map:

  an optional `list` that maps from cluster id to a cluster label, e.g.
  (1 -\> "FFA", 2 -\> "PPA")

- label:

  an optional `character` string used to label of the volume

## Value

`ClusteredNeuroVol` instance

## Details

The ClusteredNeuroVol class is designed for efficient representation and
manipulation of brain images with distinct, non-overlapping regions or
clusters. It combines the memory efficiency of sparse representations
with additional structures for managing cluster information.

The use case of `ClusteredNeuroVol` is to store volumetric data that has
been clustered into discrete sets of voxels, each of which has an
associated id. For example, this class can be used to represent
parcellated neuroimaging volumes.

## Slots

- `mask`:

  A
  [`LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.md)
  object representing the logical mask indicating the spatial domain of
  the set of clusters.

- `clusters`:

  An integer vector representing the cluster number for each voxel in
  the mask.

- `label_map`:

  A named list where each element represents a cluster and its name.

- `cluster_map`:

  An `environment` object that maps from cluster id to the set of 1D
  spatial indices belonging to that cluster.

## Methods

This class inherits methods from the
[`SparseNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/SparseNeuroVol-class.md)
class. Additional methods specific to cluster operations may be
available.

## Usage

ClusteredNeuroVol objects are particularly useful for:

- Representing parcellated brain images

- Storing results of clustering algorithms applied to neuroimaging data

- Efficient manipulation and analysis of region-based neuroimaging data

## See also

[`SparseNeuroVol-class`](https://bbuchsbaum.github.io/neuroim2/reference/SparseNeuroVol-class.md)
for the parent sparse volume class.
[`LogicalNeuroVol-class`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.md)
for the mask representation.

## Examples

``` r

# Create a simple clustered brain volume
dim <- c(10L, 10L, 10L)
mask_data <- array(rep(c(TRUE, FALSE), 500), dim)
mask <- new("LogicalNeuroVol", .Data = mask_data,
            space = NeuroSpace(dim = dim, origin = c(0,0,0), spacing = c(1,1,1)))

clusters <- as.integer(runif(sum(mask_data)) * 5)+1
label_map <- list("Cluster1" = 1, "Cluster2" = 2, "Cluster3" = 3,
                  "Cluster4" = 4, "Cluster5" = 5)

cluster_map <- list()
for (i in 1:5) {
  cluster_map[[as.character(i)]] <- which(clusters == i)
}

clustered_vol <- ClusteredNeuroVol(
                     mask = mask,
                     clusters = clusters,
                     label_map = label_map)



# Create a simple space and volume
space <- NeuroSpace(c(16, 16, 16), spacing = c(1, 1, 1))
vol_data <- array(rnorm(16^3), dim = c(16, 16, 16))
vol <- NeuroVol(vol_data, space)

# Create a binary mask (e.g., values > 0)
mask_data <- vol_data > 0
mask_vol <- LogicalNeuroVol(mask_data, space)

# Get coordinates of masked voxels
mask_idx <- which(mask_data)
coords <- index_to_coord(mask_vol, mask_idx)

# Cluster the coordinates into 10 groups
set.seed(123)  # for reproducibility
kmeans_result <- kmeans(coords, centers = 10)

# Create the clustered volume
clustered_vol <- ClusteredNeuroVol(mask_vol, kmeans_result$cluster)

# Print information about the clusters
print(clustered_vol)
#> <ClusteredNeuroVol> [10 clusters] 
#> ── Spatial ───────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 16 x 16 x 16
#>   Spacing       : 1 x 1 x 1 mm
#>   Orientation   : RAS
#> ── Clusters ──────────────────────────────────────────────────────────────────── 
#>   Count         : 10
#>   Sizes         : min=162, median=213, max=271
#>   Labels        : Clus_1, Clus_2, Clus_3, Clus_4, Clus_5, Clus_6 ...
#>   Size          : 56.8 Kb
```
