# ROIVol

A class representing a volumetric region of interest (ROI) in a brain
image, defined by a set of coordinates and associated data values.

## Slots

- `coords`:

  A `matrix` containing the 3D coordinates of the voxels within the ROI.
  Each row represents a voxel coordinate as (x, y, z).

- `.Data`:

  A `numeric` vector containing the data values associated with each
  voxel in the ROI. The length of this vector should match the number of
  rows in the `coords` matrix.

## Validity

An object of class `ROIVol` is considered valid if: - The `coords` slot
is a matrix with 3 columns. - The `.Data` slot is a numeric vector. -
The length of the `.Data` vector is equal to the number of rows in the
`coords` matrix.
