# ROIVecWindow

A class representing a spatially windowed, vector-valued volumetric
region of interest (ROI) in a brain image.

## Slots

- `coords`:

  A `matrix` containing the 3D coordinates of the voxels within the ROI.
  Each row represents a voxel coordinate as (x, y, z).

- `.Data`:

  A `matrix` containing the data values associated with each voxel in
  the ROI. Each row corresponds to a unique vector value, and the number
  of rows should match the number of rows in the `coords` matrix.

- `parent_index`:

  An `integer` representing the 1D index of the center voxel in the
  parent space.

- `center_index`:

  An `integer` representing the location in the coordinate matrix of the
  center voxel in the window.

## Validity

An object of class `ROIVecWindow` is considered valid if: - The `coords`
slot is a matrix with 3 columns. - The `.Data` slot is a matrix. - The
number of rows in the `.Data` matrix is equal to the number of rows in
the `coords` matrix.
