# Select a Subset of Clusters

Return a new object containing only the requested clusters. Clusters can
be identified by integer ID or by name (matched against the label map).

## Usage

``` r
sub_clusters(x, ids, ...)

# S4 method for class 'ClusteredNeuroVec,integer'
sub_clusters(x, ids, ...)

# S4 method for class 'ClusteredNeuroVec,numeric'
sub_clusters(x, ids, ...)

# S4 method for class 'ClusteredNeuroVec,character'
sub_clusters(x, ids, ...)

# S4 method for class 'ClusteredNeuroVol,integer'
sub_clusters(x, ids, ...)

# S4 method for class 'ClusteredNeuroVol,numeric'
sub_clusters(x, ids, ...)

# S4 method for class 'ClusteredNeuroVol,character'
sub_clusters(x, ids, ...)
```

## Arguments

- x:

  A clustered neuroimaging object.

- ids:

  Integer cluster IDs, numeric (coerced to integer), or character
  cluster names to retain.

- ...:

  Additional arguments (currently unused).

## Value

An object of the same class as `x` containing only the selected
clusters.

## Examples

``` r
sp <- NeuroSpace(c(10L, 10L, 10L), c(1, 1, 1))
mask <- LogicalNeuroVol(array(c(rep(TRUE, 500), rep(FALSE, 500)),
                              c(10, 10, 10)), sp)
clusters <- rep(1:5, length.out = 500)
cvol <- ClusteredNeuroVol(mask, clusters,
          label_map = list(A = 1, B = 2, C = 3, D = 4, E = 5))
# By integer ID
sub <- sub_clusters(cvol, c(1L, 3L))
# By name
sub2 <- sub_clusters(cvol, c("A", "C"))
```
