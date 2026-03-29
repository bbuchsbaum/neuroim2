# Space utility functions

Utilities for reasoning about mapped voxel spaces and slice embeddings.

Returns a voxel space with positive diagonal affine that encloses all
mapped input voxel centers.

Returns the affine mapping slice coordinates to 3D volume coordinates.

## Usage

``` r
output_aligned_space(mapped_voxels, voxel_sizes = NULL)

vox2out_vox(mapped_voxels, voxel_sizes = NULL)

slice_to_volume_affine(index, axis, shape = NULL, index_base = c("R", "zero"))

slice2volume(index, axis, shape = NULL, index_base = c("R", "zero"))
```

## Arguments

- mapped_voxels:

  A \`NeuroSpace\`, \`NeuroVol\`, \`NeuroVec\`, object with \`shape\`
  and \`affine\`, or length-2 sequence \`(shape, affine)\`.

- voxel_sizes:

  Optional output voxel sizes for spatial axes. If scalar, treated as
  isotropic.

- index:

  Slice index.

- axis:

  Slice axis (\`1..3\` or \`0..2\`).

- shape:

  Optional volume shape for bounds validation.

- index_base:

  Either \`"R"\` (1-based, default) or \`"zero"\`.

## Value

A named list with:

- \`shape\`: output spatial shape.

- \`affine\`: output 4x4 affine (positive diagonal).

- \`bounds\`: world-space min/max of mapped corners.

Same as \`output_aligned_space()\`.

A \`4 x 3\` affine matrix from slice coordinates to volume coordinates.

Same as \`slice_to_volume_affine()\`.

## Details

Compared to NiBabel-style helpers, these functions add a few R-friendly
improvements:

- Accept \`NeuroSpace\`, \`NeuroVol\`, and list/object \`(shape,
  affine)\` inputs.

- Handle inputs with more than 3 dimensions by using first 3 spatial
  dims.

- Support both R-style (1-based) and zero-based slice indexing.

## Examples

``` r
sp <- NeuroSpace(c(10L, 8L, 6L), spacing = c(2, 2, 2))
out <- output_aligned_space(sp)
out$shape
#> [1] 19 15 11
out$affine
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    0    0    0
#> [2,]    0    1    0    0
#> [3,]    0    0    1    0
#> [4,]    0    0    0    1

slice_aff <- slice_to_volume_affine(index = 3, axis = 3, shape = c(10, 8, 6))
slice_aff
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    2
#> [4,]    0    0    1
```
