# NeuroSpace Class

The NeuroSpace class represents the geometric properties of a brain
image, including its dimensions, origin, spacing, axes, and coordinate
transformations. It provides a comprehensive framework for handling
spatial information in neuroimaging data analysis.

## Slots

- `dim`:

  An integer vector representing the grid dimensions of the image.

- `origin`:

  A numeric vector representing the coordinates of the spatial origin.

- `spacing`:

  A numeric vector representing the dimensions (in mm) of the grid units
  (voxels).

- `axes`:

  A named
  [`AxisSet`](https://bbuchsbaum.github.io/neuroim2/reference/AxisSet-class.md)
  object representing the set of spatial axes in the untransformed
  native grid space.

- `trans`:

  A matrix representing an affine transformation that converts grid
  coordinates to real-world coordinates.

- `inverse`:

  A matrix representing an inverse transformation that converts
  real-world coordinates to grid coordinates.

## Validity

A `NeuroSpace` object is considered valid if:

- The length of the `dim` slot is equal to the lengths of the `spacing`,
  `origin`, and number of axes in the `axes` slots.

- The `dim` slot contains only non-negative values.

## Methods

The following methods are available for `NeuroSpace` objects:

- [`dim`](https://rdrr.io/r/base/dim.html): Get the dimensions of the
  space.

- [`origin`](https://bbuchsbaum.github.io/neuroim2/reference/origin-methods.md):
  Get or set the origin of the space.

- [`spacing`](https://bbuchsbaum.github.io/neuroim2/reference/spacing-methods.md):
  Get or set the spacing of the space.

- [`axes`](https://bbuchsbaum.github.io/neuroim2/reference/axes-methods.md):
  Get the axes of the space.

- [`trans`](https://bbuchsbaum.github.io/neuroim2/reference/trans-methods.md):
  Apply the affine transformation to coordinates.

## Usage

The `NeuroSpace` class is fundamental in representing and manipulating
the spatial properties of neuroimaging data. It is used extensively
throughout the package for operations that require spatial information,
such as image registration, resampling, and coordinate transformations.

## References

For more information on spatial transformations in neuroimaging: Brett,
M., Johnsrude, I. S., & Owen, A. M. (2002). The problem of functional
localization in the human brain. Nature Reviews Neuroscience, 3(3),
243-249.

## See also

[`AxisSet-class`](https://bbuchsbaum.github.io/neuroim2/reference/AxisSet-class.md)
for details on the axis set representation.
[`NeuroVol-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
and
[`NeuroVec-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md)
for classes that use NeuroSpace.

## Examples

``` r
# Create a NeuroSpace object
space <- NeuroSpace(dim = c(64L, 64L, 64L),
                    origin = c(0, 0, 0),
                    spacing = c(1, 1, 1))

# Get the dimensions
dim(space)
#> [1] 64 64 64


```
