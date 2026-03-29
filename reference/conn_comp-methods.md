# Connected components

Find connected components in an image. This function identifies and
labels spatially connected regions in neuroimaging data, supporting both
binary masks and thresholded volumes.

## Usage

``` r
conn_comp(x, ...)

# S4 method for class 'NeuroVol'
conn_comp(
  x,
  threshold = 0,
  cluster_table = TRUE,
  local_maxima = TRUE,
  local_maxima_dist = 15,
  ...
)
```

## Arguments

- x:

  the image object

- ...:

  additional arguments including:

  - threshold - numeric value defining lower intensity bound for image
    mask

  - cluster_table - logical indicating whether to return cluster
    statistics

  - local_maxima - logical indicating whether to compute local maxima

  - local_maxima_dist - minimum distance between local maxima

  - connect - connectivity pattern ("26-connect", "18-connect", or
    "6-connect")

- threshold:

  threshold defining lower intensity bound for image mask

- cluster_table:

  return cluster_table

- local_maxima:

  return table of local maxima

- local_maxima_dist:

  the distance used to define minum distance between local maxima

## Value

A list containing:

- index - A `ClusteredNeuroVol` object with cluster labels

- size - A `NeuroVol` object with cluster sizes

- voxels - A list of cluster voxel coordinates

- cluster_table - (optional) Data frame with cluster statistics

- local_maxima - (optional) Matrix of local maxima coordinates

An object representing the connected components of `x`.

## Examples

``` r
# Create a simple 3D volume with two distinct regions
space <- NeuroSpace(c(10,10,10), c(1,1,1))
vol_data <- array(0, c(10,10,10))

# Create first cluster in corner (2x2x2)
vol_data[1:2, 1:2, 1:2] <- 1

# Create second cluster in opposite corner (2x2x2)
vol_data[8:9, 8:9, 8:9] <- 1

# Create NeuroVol object
vol <- NeuroVol(vol_data, space)

# Find connected components with default 26-connectivity
# Returns components above threshold 0
comps <- conn_comp(vol, threshold=0)

# Access results
max(comps$index) == 2  # Should have 2 clusters
#> [1] TRUE
all(comps$size >= 0)    # All clusters should have >= 0
#> [1] TRUE

# Get cluster statistics
comps <- conn_comp(vol, threshold=0, cluster_table=TRUE)
# cluster_table contains: index, x, y, z, N (size), Area, value

# Find local maxima within clusters
comps <- conn_comp(vol, threshold=0, local_maxima=TRUE,
                  local_maxima_dist=2)
# local_maxima contains: index, x, y, z, value
```
