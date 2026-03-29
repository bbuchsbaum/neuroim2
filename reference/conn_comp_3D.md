# Extract Connected Components from a 3D Binary Mask

Identifies and labels connected components in a 3D binary mask using a
two-pass algorithm. The function supports different connectivity
constraints and returns both component indices and their sizes.

## Usage

``` r
conn_comp_3D(mask, connect = c("26-connect", "18-connect", "6-connect"))
```

## Arguments

- mask:

  A 3D logical array representing the binary mask

- connect:

  A character string specifying the connectivity constraint. One of
  "26-connect" (default), "18-connect", or "6-connect"

## Value

A list with the following components:

- index:

  A 3D array of integers. Each non-zero value represents the cluster
  index of the connected component for that voxel. Zero values indicate
  background.

- size:

  A 3D array of integers. Each non-zero value represents the size
  (number of voxels) of the connected component that the voxel belongs
  to. Zero values indicate background.

## Details

The function implements an efficient two-pass connected component
labeling algorithm:

- First pass: Assigns provisional labels and builds an equivalence table
  using a union-find data structure for label resolution

- Second pass: Resolves label conflicts and assigns final component
  labels

The connectivity options determine which voxels are considered adjacent:

- 6-connect: Only face-adjacent voxels (±1 step along each axis)

- 18-connect: Face and edge-adjacent voxels

- 26-connect: Face, edge, and vertex-adjacent voxels (all neighbors in a
  3x3x3 cube)

Time complexity is O(n) where n is the number of voxels in the mask,
with additional O(k) space for the union-find data structure where k is
the number of provisional labels.

## References

Rosenfeld, A., & Pfaltz, J. L. (1966). Sequential operations in digital
picture processing. Journal of the ACM, 13(4), 471-494.

## See also

[`array`](https://rdrr.io/r/base/array.html) for creating 3D arrays,
[`ClusteredNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/ClusteredNeuroVol-class.md)
for working with clustered neuroimaging data

## Examples

``` r
# Create a simple 3D binary mask with two disconnected components
mask <- array(FALSE, c(4, 4, 4))
mask[1:2, 1:2, 1:2] <- TRUE  # First component
mask[3:4, 3:4, 3:4] <- TRUE  # Second component

# Extract components using different connectivity patterns
comps <- conn_comp_3D(mask, connect = "6-connect")

# Number of components
max_comps <- max(comps$index)
cat("Found", max_comps, "components\n")
#> Found 2 components

# Size of each component
unique_sizes <- unique(comps$size[comps$size > 0])
cat("Component sizes:", paste(unique_sizes, collapse=", "), "\n")
#> Component sizes: 8 

# Try with different connectivity
comps_26 <- conn_comp_3D(mask, connect = "26-connect")
cat("Number of components with 26-connectivity:", max(comps_26$index), "\n")
#> Number of components with 26-connectivity: 1 
```
