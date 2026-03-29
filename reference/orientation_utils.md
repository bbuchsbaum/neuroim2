# Orientation utility functions

Helper functions for converting between affine matrices, axis-code
representations, and array orientation transforms.

Returns the orientation transform that maps from \`start_ornt\` to
\`end_ornt\`.

Applies flips and permutations implied by \`ornt\` to the first \`n\`
dimensions of \`arr\`, where \`n = nrow(ornt)\`.

Returns the affine mapping transformed array coordinates back to source
array coordinates for a transform encoded by \`ornt\`.

## Usage

``` r
affine_to_orientation(affine, tol = NULL)

orientation_transform(start_ornt, end_ornt)

apply_orientation(arr, ornt)

orientation_inverse_affine(ornt, shape)

orientation_to_axcodes(ornt, labels = NULL)

axcodes_to_orientation(axcodes, labels = NULL)

affine_to_axcodes(affine, labels = NULL, tol = NULL)
```

## Arguments

- affine:

  Affine matrix.

- tol:

  Optional singular-value tolerance.

- start_ornt:

  Starting orientation matrix.

- end_ornt:

  Target orientation matrix.

- arr:

  Array-like object.

- ornt:

  Orientation matrix.

- shape:

  Shape of source array.

- labels:

  Optional label pairs per axis.

- axcodes:

  Character vector of axis codes. \`NA\` indicates dropped axis.

## Value

A \`p x 2\` orientation matrix with columns \`axis\` and \`flip\`.

Orientation matrix representing \`start -\> end\`.

Reoriented array.

Homogeneous affine matrix of size \`(p + 1) x (p + 1)\`.

Character vector of axis codes (positive ends), with \`NA\` for dropped
axes.

Orientation matrix.

Character vector of axis codes.

## Details

Orientation matrices (\`ornt\`) use two columns:

- column 1 (\`axis\`) stores the output axis index (1-based),

- column 2 (\`flip\`) stores direction (\`1\` or \`-1\`).

This is an R counterpart to NiBabel's orientation utilities, adapted to
\`neuroim2\` conventions.

## Examples

``` r
aff <- diag(4)
ornt <- affine_to_orientation(aff)
orientation_to_axcodes(ornt)
#> [1] "R" "A" "S"

arr <- array(1:24, dim = c(2, 3, 4))
tx <- orientation_transform(
  axcodes_to_orientation(c("R", "A", "S")),
  axcodes_to_orientation(c("A", "R", "S"))
)
out <- apply_orientation(arr, tx)

inv_aff <- orientation_inverse_affine(tx, dim(arr))
inv_aff
#>      [,1] [,2] [,3] [,4]
#> [1,]    0    1    0    0
#> [2,]    1    0    0    0
#> [3,]    0    0    1    0
#> [4,]    0    0    0    1
```
