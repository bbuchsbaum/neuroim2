# ROIVolWindow

A class representing a spatially windowed volumetric region of interest
(ROI) in a brain image, derived from a larger parent ROI.

## Slots

- `parent_index`:

  An `integer` representing the 1D index of the center voxel in the
  parent space.

- `center_index`:

  An `integer` representing the location in the coordinate matrix of the
  center voxel in the window.

- `coords`:

  A `matrix` containing the 3D coordinates of the voxels within the ROI.
  Each row represents a voxel coordinate as (x, y, z).

- `.Data`:

  A `numeric` vector containing the data values associated with each
  voxel in the ROI. The length of this vector should match the number of
  rows in the `coords` matrix.

## Validity

An object of class `ROIVolWindow` is considered valid if: - The `coords`
slot is a matrix with 3 columns. - The `.Data` slot is a numeric
vector. - The length of the `.Data` vector is equal to the number of
rows in the `coords` matrix.
