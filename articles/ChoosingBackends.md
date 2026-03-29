# Choosing a Backend

## Why Backends Matter

`neuroim2` supports multiple representations for the same conceptual
data: 3D volumes, 4D time-series, sparse masked data, memory-mapped
arrays, and file-backed access. Choosing the right representation is
usually the biggest practical decision you make up front.

This guide is the shortest path to that choice.

## Quick Rules

- Use `DenseNeuroVol` or `DenseNeuroVec` when the data fits comfortably
  in memory and you want the simplest behavior.
- Use `SparseNeuroVol` or `SparseNeuroVec` when only a subset of voxels
  is meaningful and absent voxels should be treated as missing support
  rather than literal zeros.
- Use `MappedNeuroVec` when you want fast random access to an
  uncompressed NIfTI on disk without loading it all into memory.
- Use `FileBackedNeuroVec` when you want on-demand volume access from a
  file and can tolerate more I/O per access.
- Use `NeuroHyperVec` when the input is genuinely 5D.

## Decision Table

| Situation                                         | Recommended backend  | Why                             |
|:--------------------------------------------------|:---------------------|:--------------------------------|
| Small-to-moderate 3D image                        | `DenseNeuroVol`      | simplest, lowest overhead       |
| Small-to-moderate 4D fMRI series                  | `DenseNeuroVec`      | easiest indexing and arithmetic |
| Masked analysis over a subset of voxels           | `SparseNeuroVec`     | stores only supported voxels    |
| Large uncompressed NIfTI with random-access needs | `MappedNeuroVec`     | memory mapping is efficient     |
| Large file, on-demand volume access               | `FileBackedNeuroVec` | no full materialization         |
| 5D image                                          | `NeuroHyperVec`      | intended representation         |

## Dense In-Memory Data

If the dataset fits in RAM, start here.

``` r
file_name <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
vec <- read_vec(file_name)

class(vec)
#> [1] "DenseNeuroVec"
#> attr(,"package")
#> [1] "neuroim2"
dim(vec)
#> [1] 64 64 25  4
```

This gives you the most predictable behavior for:

- arithmetic
- indexing
- resampling
- plotting
- converting to matrices or arrays

## Sparse Masked Data

If you only care about a subset of voxels, read the data through a mask
and keep the representation sparse.

``` r
file_name <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
mask_vol <- read_vol(file_name) > 0
svec <- read_vec(file_name, mask = mask_vol)

class(svec)
#> [1] "SparseNeuroVec"
#> attr(,"package")
#> [1] "neuroim2"
dim(svec)
#> [1] 64 64 25  4
sum(mask(svec))
#> [1] 29532
```

Use this when:

- your workflow is defined on a brain mask or ROI support
- you want to avoid materializing out-of-mask voxels
- absent voxels should be treated as missing support

## Memory-Mapped Access

For large uncompressed NIfTI files, `mmap` mode gives efficient on-disk
access.

``` r
file_name <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
mvec <- read_vec(file_name, mode = "mmap")

series(mvec, 1, 1, 1)
sub_vector(mvec, 1:5)
```

Use this when:

- the file is uncompressed
- you need repeated random access
- you want to avoid holding the full array in memory

## File-Backed Access

`filebacked` mode is the safest on-disk choice when you want deferred
access without requiring the full dense object in memory.

``` r
file_name <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
fbvec <- read_vec(file_name, mode = "filebacked")

fbvec[[1]]
sub_vector(fbvec, 1:3)
```

Use this when:

- you mostly work volume-by-volume
- you want an explicit on-disk representation
- you do not want to materialize the entire dataset

## 5D Data

When the input is 5D, use
[`read_image()`](https://bbuchsbaum.github.io/neuroim2/reference/read_image.md)
or
[`read_hyper_vec()`](https://bbuchsbaum.github.io/neuroim2/reference/read_hyper_vec.md).

``` r
img5d <- read_image("some_5d_image.nii.gz")
class(img5d)
```

## Recommended Defaults

If you are unsure:

1.  Start with
    [`read_vec()`](https://bbuchsbaum.github.io/neuroim2/reference/read_vec.md)
    or
    [`read_vol()`](https://bbuchsbaum.github.io/neuroim2/reference/read_vol.md).
2.  Move to `SparseNeuroVec` only when a mask meaningfully reduces the
    support.
3.  Move to `mmap` or `filebacked` when memory or I/O costs force the
    issue.

## Related Articles

- [`vignette("Overview")`](https://bbuchsbaum.github.io/neuroim2/articles/Overview.md)
  for general orientation
- [`vignette("NeuroVector")`](https://bbuchsbaum.github.io/neuroim2/articles/NeuroVector.md)
  for 4D workflows
- [`vignette("Resampling")`](https://bbuchsbaum.github.io/neuroim2/articles/Resampling.md)
  for spatial transforms
- [`vignette("regionOfInterest")`](https://bbuchsbaum.github.io/neuroim2/articles/regionOfInterest.md)
  for ROI and searchlight workflows
