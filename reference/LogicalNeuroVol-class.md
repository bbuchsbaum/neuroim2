# LogicalNeuroVol Class

This class represents a three-dimensional brain image where all values
are either TRUE or FALSE. It is particularly useful for creating and
managing binary masks for brain images.

This function constructs a `LogicalNeuroVol` instance.

## Usage

``` r
LogicalNeuroVol(data, space, label = "", indices = NULL)
```

## Arguments

- data:

  A three-dimensional `array`, a 1D vector with length equal to
  `prod(dim(space))`, or a set of `indices` where elements are `TRUE`.

- space:

  An instance of class
  [`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md).

- label:

  A `character` string.

- indices:

  An optional 1-d index vector.

## Value

A `LogicalNeuroVol` instance.

## Details

The LogicalNeuroVol class extends the
[`DenseNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/DenseNeuroVol-class.md)
class, inheriting its spatial properties and array-based storage.
However, it constrains the values to be logical (TRUE or FALSE), making
it ideal for representing binary masks, regions of interest (ROIs), or
segmentation results in neuroimaging analyses.

## Slots

- `.Data`:

  A logical array containing the binary volume data.

- `space`:

  A
  [`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md)
  object defining the spatial properties of the volume.

## Methods

This class inherits methods from
[`DenseNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/DenseNeuroVol-class.md).
Additional methods specific to logical operations may be available.

## See also

[`DenseNeuroVol-class`](https://bbuchsbaum.github.io/neuroim2/reference/DenseNeuroVol-class.md)
for the parent class.
[`NeuroVol-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
for the base volumetric image class.

## Examples

``` r
# Create a simple logical brain volume (e.g., a mask)
dim <- c(64L, 64L, 64L)
mask_data <- array(sample(c(TRUE, FALSE), prod(dim), replace = TRUE), dim)
mask_space <- NeuroSpace(dim = dim, origin = c(0, 0, 0), spacing = c(1, 1, 1))
brain_mask <- new("LogicalNeuroVol", .Data = mask_data, space = mask_space)

# Check the proportion of TRUE voxels
true_proportion <- sum(brain_mask) / prod(dim(brain_mask))
print(paste("Proportion of TRUE voxels:", true_proportion))
#> [1] "Proportion of TRUE voxels: 0.500423431396484"

# Load an example brain mask
brain_mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))

# Convert the brain mask to a LogicalNeuroVol
logical_vol <- LogicalNeuroVol(brain_mask, space(brain_mask))
```
