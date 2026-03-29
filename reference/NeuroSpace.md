# NeuroSpace: Spatial Reference System for Neuroimaging Data

The `NeuroSpace` class defines the spatial properties and coordinate
system of neuroimaging data. It encapsulates all information needed to
map between voxel indices and real-world coordinates, including
dimensions, voxel spacing, origin, axis orientation, and coordinate
transformations.

## Usage

``` r
NeuroSpace(dim, spacing = NULL, origin = NULL, axes = NULL, trans = NULL)
```

## Arguments

- dim:

  An integer vector specifying the dimensions of the image grid. Must be
  positive.

- spacing:

  A numeric vector specifying the physical size of each voxel (typically
  in millimeters). Must be positive. If NULL, defaults to ones.

- origin:

  A numeric vector specifying the real-world coordinates of the first
  voxel. If NULL, defaults to zeros.

- axes:

  An
  [`AxisSet`](https://bbuchsbaum.github.io/neuroim2/reference/AxisSet-class.md)
  object defining the orientation and ordering of the coordinate axes.
  If NULL, defaults to standard neurological convention
  (Left-Posterior-Inferior for 3D).

- trans:

  A transformation matrix mapping voxel indices to world coordinates. If
  NULL, constructed from spacing and origin.

## Value

A new
[`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md)
object

## Details

Spatial Reference System for Neuroimaging Data

## Coordinate Systems

NeuroSpace manages two coordinate systems:

- Voxel coordinates: Zero-based indices into the image grid

- World coordinates: Real-world coordinates (typically in millimeters)

The transformation between these systems is defined by:

- Voxel spacing (physical size of voxels)

- Origin (world coordinates of first voxel)

- Axis orientation (how image axes map to anatomical directions)

## Validation

The constructor performs extensive validation:

- All dimensions must be positive integers

- All spacing values must be positive

- Origin and spacing must have matching lengths

- Transformation matrix must be invertible

## References

For details on neuroimaging coordinate systems:

- Brett, M., Johnsrude, I. S., & Owen, A. M. (2002). The problem of
  functional localization in the human brain. Nature Reviews
  Neuroscience, 3(3), 243-249.

- Evans, A. C., et al. (1993). 3D statistical neuroanatomical models
  from 305 MRI volumes. Nuclear Science Symposium and Medical Imaging
  Conference.

## See also

[`AxisSet`](https://bbuchsbaum.github.io/neuroim2/reference/AxisSet-class.md)
for axis orientation specification,
[`coord_to_index`](https://bbuchsbaum.github.io/neuroim2/reference/coord_to_index-methods.md)
for coordinate conversion,
[`index_to_coord`](https://bbuchsbaum.github.io/neuroim2/reference/index_to_coord-methods.md)
for inverse coordinate conversion,
[`NeuroObj`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroObj-class.md)
for objects using NeuroSpace

## Examples

``` r
# Create a standard 3D space (64x64x40 voxels, 2mm isotropic)
space_3d <- NeuroSpace(
  dim = c(64L, 64L, 40L),
  spacing = c(2, 2, 2),
  origin = c(-90, -126, -72)
)

# Check properties
dim(space_3d)           # Image dimensions
#> [1] 64 64 40
spacing(space_3d)       # Voxel sizes
#> [1] 2 2 2
origin(space_3d)        # World-space origin
#> [1]  -90 -126  -72

# Create a 2D slice space
space_2d <- NeuroSpace(
  dim = c(128L, 128L),
  spacing = c(1.5, 1.5),
  origin = c(-96, -96)
)

# Convert between coordinate systems
world_coords <- c(0, 0, 0)
vox_idx <- coord_to_index(space_3d, world_coords)
back_to_world <- index_to_coord(space_3d, vox_idx)
```
