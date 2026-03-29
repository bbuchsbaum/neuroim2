# Constructor for NeuroHyperVec class

Constructor for NeuroHyperVec class

## Usage

``` r
NeuroHyperVec(data, space, mask)
```

## Arguments

- data:

  A matrix or three-dimensional array containing the data.

- space:

  A
  [`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md)
  object defining the spatial dimensions.

- mask:

  A mask volume (array, vector, or
  [`LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.md)).

## Value

A new
[`NeuroHyperVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroHyperVec-class.md)
object.

## See also

[`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md),
[`LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.md)

## Examples

``` r
# Create a 5D space (10x10x10 spatial, 2 trials, 2 features)
space <- NeuroSpace(c(10,10,10,2,2))

# Create a mask for the spatial dimensions
space3d <- NeuroSpace(c(10,10,10))
mask_data <- array(TRUE, dim=c(10,10,10))  # All voxels active
mask <- LogicalNeuroVol(mask_data, space3d)

# Create data in the format expected by NeuroHyperVec:
# 3D array with dimensions [features x trials x voxels]
n_features <- 2
n_trials <- 2
n_voxels <- sum(mask_data)  # 1000 voxels
data_array <- array(rnorm(n_features * n_trials * n_voxels),
                   dim = c(n_features, n_trials, n_voxels))

# Create the NeuroHyperVec object
hvec <- NeuroHyperVec(data_array, space, mask)
```
