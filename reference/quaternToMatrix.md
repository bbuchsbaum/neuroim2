# Convert Quaternion Parameters to a Transformation Matrix

Given a quaternion `(b, c, d)`, a scalar offset (origin), voxel step
sizes, and the `qfac` sign, reconstructs a 4x4 affine matrix
representing rotation, scaling, and translation as used in NIfTI-1.

## Usage

``` r
quaternToMatrix(quat, origin, stepSize, qfac)
```

## Arguments

- quat:

  A numeric vector of length 3 containing the quaternion parameters
  \\(b, c, d)\\. The scalar part \\a\\ is computed internally.

- origin:

  A numeric vector of length 3 specifying the translation components
  (often the real-space origin or offset).

- stepSize:

  A numeric vector of length 3 giving the voxel dimensions along each
  axis (e.g., `(dx, dy, dz)`).

- qfac:

  Either `+1` or `-1`, indicating the sign from the determinant check in
  [`matrixToQuatern`](https://bbuchsbaum.github.io/neuroim2/reference/matrixToQuatern.md).

## Value

A 4x4 numeric affine transformation matrix. The top-left 3x3 submatrix
encodes rotation and scaling, and the 4th column encodes translation.

## Details

This function uses the quaternion formalism common in neuroimaging,
adding the offset (translation) into the 4th column, and applying the
voxel sizes along each axis. If `qfac` is `-1`, the \\z\\ scale is
negated. The resulting 4x4 matrix is typically used as an affine
transform for voxel-to-world coordinate mapping.

## See also

[`matrixToQuatern`](https://bbuchsbaum.github.io/neuroim2/reference/matrixToQuatern.md)
for converting a matrix back to quaternion form.
