# NeuroSlice: 2D Neuroimaging Data Container

Creates a `NeuroSlice` object representing a two-dimensional slice of
neuroimaging data with associated spatial information. This class is
particularly useful for working with individual slices from volumetric
neuroimaging data or for visualizing 2D cross-sections.

## Usage

``` r
NeuroSlice(data, space, indices = NULL)
```

## Arguments

- data:

  A vector or matrix containing the slice data values.

- space:

  An object of class
  [`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md)
  defining the spatial properties (dimensions, spacing, origin) of the
  slice.

- indices:

  Optional integer vector. When `data` is provided as a 1D vector,
  `indices` specifies the linear indices where the data values should be
  placed in the 2D slice. Useful for creating sparse slices. Default is
  NULL.

## Value

A new object of class
[`NeuroSlice`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSlice-class.md).

## Details

Two-Dimensional Neuroimaging Data Slice

## Input Validation

The function performs several validation checks:

- Verifies that `space` is 2-dimensional

- Ensures data dimensions are compatible with `space`

- Validates `indices` when provided for sparse initialization

## Data Handling

The function supports two initialization modes:

- Dense mode (indices = NULL):

  - Data is reshaped if necessary to match space dimensions

  - Dimensions must match exactly after reshaping

- Sparse mode (indices provided):

  - Creates a zero-initialized matrix matching space dimensions

  - Places data values at specified indices

## See also

[`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md)
for defining spatial properties,
[`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
for 3D volumetric data,
[`plot`](https://rdrr.io/r/graphics/plot.default.html) for visualization
methods

## Examples

``` r
# Create a 64x64 slice space
slice_space <- NeuroSpace(c(64, 64), spacing = c(2, 2))

# Example 1: Dense slice from matrix
slice_data <- matrix(rnorm(64*64), 64, 64)
dense_slice <- NeuroSlice(slice_data, slice_space)

# Example 2: Dense slice from vector
vec_data <- rnorm(64*64)
vec_slice <- NeuroSlice(vec_data, slice_space)

# Example 3: Sparse slice with specific values
n_points <- 100
sparse_data <- rnorm(n_points)
sparse_indices <- sample(1:(64*64), n_points)
sparse_slice <- NeuroSlice(sparse_data, slice_space, indices = sparse_indices)
```
