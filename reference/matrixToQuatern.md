# Convert a Transformation Matrix to a Quaternion Representation

Extracts the rotation and scaling components from a 3x3 (or 4x4)
transformation matrix, normalizes them, and computes the corresponding
quaternion parameters and a sign factor (\`qfac\`) indicating whether
the determinant is negative.

## Usage

``` r
matrixToQuatern(mat)
```

## Arguments

- mat:

  A numeric matrix with at least the top-left 3x3 portion containing
  rotation/scaling. Often a 4x4 affine transform, but only the 3x3
  top-left submatrix is used in practice.

## Value

A named `list` with two elements:

- `quaternion`:

  A numeric vector of length 3, \\(b, c, d)\\, which—together with \\a\\
  derived internally—represents the rotation.

- `qfac`:

  Either `+1` or `-1`, indicating whether the determinant of the
  rotation submatrix is positive or negative, respectively.

## Details

This function first checks and corrects for zero-length axes in the
upper-left corner of the matrix, then normalizes each column to extract
the pure rotation. If the determinant of the rotation submatrix is
negative, the `qfac` is set to `-1`, and the third column is negated.
Finally, the quaternion parameters (\\a, b, c, d\\) are computed
following standard NIfTI-1 conventions for representing the rotation in
3D.

## References

\- Cox RW. \*Analysis of Functional NeuroImages\* (AFNI) and NIfTI-1
quaternion conventions. <https://afni.nimh.nih.gov>

## See also

[`quaternToMatrix`](https://bbuchsbaum.github.io/neuroim2/reference/quaternToMatrix.md)
for the inverse operation, converting quaternion parameters back to a
transform matrix.
