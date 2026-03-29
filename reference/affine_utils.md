# Affine utility functions

Utilities for point transforms and affine decomposition/composition.

Applies a homogeneous affine transform to points where the coordinate
axis is the last dimension.

Useful when expanding an affine to include additional dimensions.

For arrays \`A, B, C, ...\`, returns \`A

Returns per-output-axis obliquity relative to cardinal axes, in radians.

Preserves rotations/shears and world location of the center voxel.

## Usage

``` r
apply_affine(aff, pts, inplace = FALSE)

to_matvec(transform)

from_matvec(matrix, vector = NULL)

append_diag(aff, steps, starts = numeric(0))

dot_reduce(...)

voxel_sizes(affine)

obliquity(affine)

rescale_affine(affine, shape, zooms, new_shape = NULL)
```

## Arguments

- aff:

  Base affine matrix.

- pts:

  Points with last dimension \`N - 1\`; shape may be vector, matrix, or
  higher-dimensional array.

- inplace:

  Included for API compatibility; ignored in R.

- transform:

  Homogeneous transform matrix.

- matrix:

  Linear part of transform (\`N x M\`).

- vector:

  Optional translation vector of length \`N\`.

- steps:

  Diagonal values for appended dimensions.

- starts:

  Optional translations for appended dimensions.

- ...:

  Matrices/vectors compatible with `%*%`.

- affine:

  Square homogeneous affine matrix.

- shape:

  Original grid shape (spatial axes).

- zooms:

  New voxel sizes for spatial axes.

- new_shape:

  Optional new grid shape; defaults to \`shape\`.

## Value

Points transformed by \`aff\`, with same leading shape as \`pts\`.

A list with \`matrix\` and \`vector\`.

Homogeneous affine matrix of shape \`(N+1) x (M+1)\`.

Expanded affine matrix.

Matrix product.

Numeric vector of voxel sizes (column norms of linear block).

Numeric vector of obliquity angles.

Rescaled affine matrix.

## Details

These functions mirror core NiBabel affine helpers while using R-first
conventions and stronger argument checks.

## Examples

``` r
aff <- diag(c(2, 3, 4, 1))
aff[1:3, 4] <- c(10, 20, 30)

pts <- rbind(c(1, 2, 3), c(4, 5, 6))
apply_affine(aff, pts)
#>      [,1] [,2] [,3]
#> [1,]   12   26   42
#> [2,]   18   35   54

mv <- to_matvec(aff)
from_matvec(mv$matrix, mv$vector)
#>      [,1] [,2] [,3] [,4]
#> [1,]    2    0    0   10
#> [2,]    0    3    0   20
#> [3,]    0    0    4   30
#> [4,]    0    0    0    1
```
