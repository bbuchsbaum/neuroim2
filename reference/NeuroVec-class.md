# NeuroVec Class

This S4 class represents a four-dimensional brain image, which is used
to store and process time series neuroimaging data such as fMRI or 4D
functional connectivity maps. The class extends the basic functionality
of
[`NeuroObj`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroObj-class.md).

The `NeuroVec` class represents a vectorized form of neuroimaging data,
supporting both in-memory and file-backed data modes. It provides
efficient data storage and access methods and integrates with the
spatial reference system provided by
[`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md).

## Usage

``` r
NeuroVec(data, space = NULL, mask = NULL, label = "")
```

## Arguments

- data:

  The image data. This can be:

  - A matrix (voxels x time points)

  - A 4D array

  - A list of
    [`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
    objects

  If a list of NeuroVol objects is provided, the geometric space
  (`NeuroSpace`) will be inferred from the constituent volumes, which
  must all be identical.

- space:

  An optional
  [`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md)
  object defining the spatial properties of the image. Not required if
  `data` is a list of NeuroVol objects.

- mask:

  An optional logical array specifying which voxels to include. If
  provided, a SparseNeuroVec object will be created.

- label:

  A character string providing a label for the NeuroVec object. Default
  is an empty string.

## Value

A concrete instance of the `NeuroVec` class:

- If `mask` is provided: a
  [`SparseNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/SparseNeuroVec-class.md)
  object

- Otherwise: a
  [`DenseNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/DenseNeuroVec-class.md)
  object

## Details

NeuroVec objects are designed to handle 4D neuroimaging data, where the
first three dimensions represent spatial coordinates, and the fourth
dimension typically represents time or another series dimension. This
structure is particularly useful for storing and analyzing functional
MRI data, time series of brain states, or multiple 3D volumes in a
single object.

The function performs several operations:

- If `data` is a list of NeuroVol objects, it combines them into a
  single 4D array.

- It checks that the dimensions of `data` match the provided `space`.

- Depending on whether a `mask` is provided, it creates either a
  DenseNeuroVec or a SparseNeuroVec object.

## Slots

- space:

  A
  [`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md)
  object defining the spatial properties of the image.

- label:

  A character string providing a label for the NeuroVec object.

## Methods

Methods specific to NeuroVec objects may include operations for time
series analysis, 4D data manipulation, and extraction of 3D volumes or
time courses.

## Usage

To create a NeuroVec object, use the constructor function `NeuroVec()`.
This function should handle the appropriate initialization of the 4D
data structure and associated spatial information.

## See also

[`NeuroObj-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroObj-class.md)
for the parent class.
[`DenseNeuroVec-class`](https://bbuchsbaum.github.io/neuroim2/reference/DenseNeuroVec-class.md)
and
[`SparseNeuroVec-class`](https://bbuchsbaum.github.io/neuroim2/reference/SparseNeuroVec-class.md)
for specific implementations.

[`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md)
for spatial information,
[`sub_vector`](https://bbuchsbaum.github.io/neuroim2/reference/sub_vector-methods.md)
for subsetting routines, and
[`index_to_coord`](https://bbuchsbaum.github.io/neuroim2/reference/index_to_coord-methods.md)
for coordinate conversion.
[`DenseNeuroVec-class`](https://bbuchsbaum.github.io/neuroim2/reference/DenseNeuroVec-class.md),
[`SparseNeuroVec-class`](https://bbuchsbaum.github.io/neuroim2/reference/SparseNeuroVec-class.md)
for the specific NeuroVec types.
[`NeuroVol-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
for 3D volumetric data.

## Examples

``` r
# Load an example 4D brain image
example_4d_image <- read_vec(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))

# Create a NeuroVec object
neuro_vec <- NeuroVec(data = array(rnorm(64*64*32*10), dim = c(64, 64, 32, 10)),
                      space = NeuroSpace(dim = c(64, 64, 32,10),
                      origin = c(0, 0, 0),
                      spacing = c(3, 3, 4)))


dim(neuro_vec)
#> [1] 64 64 32 10

# Extract a single 3D volume (e.g., the first time point)
first_volume <- neuro_vec[[1]]


# Load an example 4D brain image
example_file <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
example_4d_image <- read_vec(example_file)

# Create a DenseNeuroVec object
dense_vec <- NeuroVec(data = example_4d_image@.Data,
                      space = space(example_4d_image))
print(dense_vec)
#> <DenseNeuroVec> [3.1 Mb] 
#> ── Spatial ───────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 64 x 64 x 25 (4 timepoints)
#>   Spacing       : 3.5 x 3.5 x 3.7
#>   Origin        : 112, -108, -46.2
#>   Orientation   : LAS
#> ── Data ──────────────────────────────────────────────────────────────────────── 
#>   Mean +/- SD   : 0.288 +/- 0.453 (t=1)

# Create a SparseNeuroVec object with a mask
mask <- array(runif(prod(dim(example_4d_image)[1:3])) > 0.5,
              dim = dim(example_4d_image)[1:3])
sparse_vec <- NeuroVec(data = example_4d_image@.Data,
                       space = space(example_4d_image),
                       mask = mask)
print(sparse_vec)
#> <SparseNeuroVec> [2.6 Mb] 
#> ── Spatial ───────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 64 x 64 x 25
#>   Time Points   : 4
#>   Spacing       : 3.5 x 3.5 x 3.7
#>   Origin        : 112, -108, -46.2
#> ── Sparse ────────────────────────────────────────────────────────────────────── 
#>   Cardinality   : 51129
```
