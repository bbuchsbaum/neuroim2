# Deoblique a Neuroimaging Space or Volume

AFNI-like helper that mirrors the core behavior of `3dWarp -deoblique`:

- If `gridset` is supplied, use that as the output grid.

- Otherwise, build an axis-aligned output grid that encloses the input
  field-of-view.

- If `newgrid` is not supplied, use the minimum input voxel size
  isotropically (AFNI-style default for deobliquing).

## Usage

``` r
deoblique(
  x,
  gridset = NULL,
  newgrid = NULL,
  method = c("linear", "nearest", "cubic"),
  engine = c("internal")
)
```

## Arguments

- x:

  A `NeuroSpace` or `NeuroVol`.

- gridset:

  Optional output grid (a `NeuroSpace` or `NeuroVol`), analogous to
  AFNI's `-gridset`. Mutually exclusive with `newgrid`.

- newgrid:

  Optional scalar output voxel size (mm), analogous to AFNI's
  `-newgrid`. If omitted and `gridset` is `NULL`, the minimum input
  voxel size is used isotropically.

- method:

  Interpolation method used when `x` is a `NeuroVol`: one of
  `"nearest"`, `"linear"`, or `"cubic"`.

- engine:

  Resampling engine passed to
  [`resample_to`](https://bbuchsbaum.github.io/neuroim2/reference/resample_to.md).

## Value

If `x` is a `NeuroSpace`, returns a deobliqued `NeuroSpace`. If `x` is a
`NeuroVol`, returns a resampled `NeuroVol` in deobliqued space.

## Details

For `NeuroSpace`, this returns the target deobliqued space. For
`NeuroVol`, it also resamples image data into that space.

## See also

[`output_aligned_space`](https://bbuchsbaum.github.io/neuroim2/reference/space_utils.md),
[`resample_to`](https://bbuchsbaum.github.io/neuroim2/reference/resample_to.md)

## Examples

``` r
sp <- NeuroSpace(c(32, 32, 20), spacing = c(2, 2, 3))
tx <- trans(sp)
tx[1, 2] <- 0.15
sp_obl <- NeuroSpace(dim(sp), spacing = spacing(sp), trans = tx)

# Build deobliqued target space (minimum spacing default)
sp_deob <- deoblique(sp_obl)

# Resample a volume to deobliqued space
vol <- NeuroVol(array(rnorm(prod(dim(sp_obl))), dim(sp_obl)), sp_obl)
# \donttest{
vol_deob <- deoblique(vol, method = "linear")
# }
```
