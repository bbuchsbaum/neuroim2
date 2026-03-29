# DenseNeuroVec Class

A class representing a four-dimensional brain image, backed by a dense
array. This class is designed for neuroimaging data where most voxels
contain non-zero values.

This function constructs a DenseNeuroVec object, which represents a
dense four-dimensional brain image. It handles various input data
formats and ensures proper dimensionality.

## Usage

``` r
DenseNeuroVec(data, space, label = "none")
```

## Arguments

- data:

  The image data. This can be:

  - A 4-dimensional array

  - A 2-dimensional matrix (either nvoxels x ntime-points or
    ntime-points x nvoxels)

  - A vector (which will be reshaped to match the space dimensions)

- space:

  A
  [`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md)
  object defining the spatial properties of the image.

- label:

  A character string providing a label for the DenseNeuroVec object.
  Default is an empty string.

## Value

A concrete instance of the `DenseNeuroVec` class.

## Details

DenseNeuroVec objects store their data in a dense array format, which is
efficient for operations that require frequent access to all voxels.
This class inherits from both
[`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md)
and `array` classes, combining spatial information with array-based
storage.

The function performs several operations based on the input data type:

- For matrix input: It determines the correct orientation (voxels x time
  or time x voxels) and reshapes accordingly. If necessary, it adds a
  4th dimension to the space object.

- For vector input: It reshapes the data to match the dimensions
  specified in the space object.

- For array input: It ensures the dimensions match those specified in
  the space object.

Note that the label parameter is currently not used in the object
creation, but is included for potential future use or consistency with
other constructors.

## Validity

A DenseNeuroVec object is considered valid if:

- The underlying data is a four-dimensional array.

## See also

[`NeuroVec-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md)
for the parent class.
[`SparseNeuroVec-class`](https://bbuchsbaum.github.io/neuroim2/reference/SparseNeuroVec-class.md)
for a sparse representation alternative.

[`NeuroVec-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md)
for the parent class.
[`SparseNeuroVec-class`](https://bbuchsbaum.github.io/neuroim2/reference/SparseNeuroVec-class.md)
for the sparse version of 4D brain images.
[`NeuroSpace-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md)
for details on spatial properties.

## Examples

``` r
# Create a simple 4D brain image
data <- array(rnorm(64*64*32*10), dim = c(64, 64, 32, 10))
space <- NeuroSpace(dim = c(64, 64, 32,10), origin = c(0, 0, 0), spacing = c(3, 3, 4))
dense_vec <- new("DenseNeuroVec", .Data = data, space = space)

# Access dimensions
dim(dense_vec)
#> [1] 64 64 32 10

# Extract a single 3D volume
first_volume <- dense_vec[[1]]


# Create a simple 4D brain image
dim <- c(64, 64, 32, 10)  # 64x64x32 volume with 10 time points
data <- array(rnorm(prod(dim)), dim)
space <- NeuroSpace(dim, spacing = c(3, 3, 4))

# Create a DenseNeuroVec object
dense_vec <- DenseNeuroVec(data = data, space = space, label = "Example")
print(dense_vec)
#> <DenseNeuroVec> [10 Mb] 
#> ── Spatial ───────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 64 x 64 x 32 (10 timepoints)
#>   Spacing       : 3 x 3 x 4
#>   Origin        : 0, 0, 0
#>   Orientation   : RAS
#> ── Data ──────────────────────────────────────────────────────────────────────── 
#>   Mean +/- SD   : 0.005 +/- 1.002 (t=1)
#>   Label         : Example

# Create from a matrix (voxels x time)
mat_data <- matrix(rnorm(prod(dim)), nrow = prod(dim[1:3]), ncol = dim[4])
dense_vec_mat <- DenseNeuroVec(data = mat_data, space = space)
print(dense_vec_mat)
#> <DenseNeuroVec> [10 Mb] 
#> ── Spatial ───────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 64 x 64 x 32 (10 timepoints)
#>   Spacing       : 3 x 3 x 4
#>   Origin        : 0, 0, 0
#>   Orientation   : RAS
#> ── Data ──────────────────────────────────────────────────────────────────────── 
#>   Mean +/- SD   : 0.002 +/- 0.999 (t=1)
#>   Label         : none
```
