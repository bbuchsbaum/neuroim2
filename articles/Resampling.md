# Resampling, downsampling, and reorientation

This vignette shows how to:

- Resample an image into a new grid
  ([`resample_to()`](https://bbuchsbaum.github.io/neuroim2/reference/resample_to.md)
  /
  [`resample()`](https://bbuchsbaum.github.io/neuroim2/reference/resample-methods.md)).
- Downsample 3D/4D data to coarser resolution
  ([`downsample()`](https://bbuchsbaum.github.io/neuroim2/reference/downsample-methods.md)).
- Change orientation while keeping physical coordinates consistent
  ([`reorient()`](https://bbuchsbaum.github.io/neuroim2/reference/reorient-methods.md)).

We’ll use the example NIfTI that ships with the package:

``` r
demo_path <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
vol  <- read_vol(demo_path)  # 3D NeuroVol
vec4 <- read_vec(demo_path)  # 4D NeuroVec

sp   <- space(vol)
dim(vol)
#> [1] 64 64 25
spacing(vol)
#> [1] 3.5 3.5 3.7
```

## 1. Resampling to a new space

[`resample_to()`](https://bbuchsbaum.github.io/neuroim2/reference/resample_to.md)
is a convenience wrapper around the S4
[`resample()`](https://bbuchsbaum.github.io/neuroim2/reference/resample-methods.md)
methods that accepts human‑readable interpolation names.

You specify a *target* grid (`NeuroSpace` or `NeuroVol`) and choose
`"nearest"`, `"linear"`, or `"cubic"` interpolation:

``` r
# Create a target space with 2× smaller voxels in each dimension
sp_fine <- NeuroSpace(
  dim    = sp@dim * 2L,
  spacing = sp@spacing / 2,
  origin  = sp@origin,
  trans   = trans(vol)
)

vol_fine_lin <- resample_to(vol, sp_fine, method = "linear")
dim(vol_fine_lin)
#> [1] 128 128  50
spacing(vol_fine_lin)
#> [1] 3.5 3.5 3.7
```

Typical patterns:

- Use `"nearest"` when resampling label or mask images.
- Use `"linear"` or `"cubic"` for continuous data (e.g. intensities,
  statistics).

You can also call `resample(source, target, interpolation = 0/1/3)`
directly if you prefer numeric codes.

## 2. Downsampling to coarser resolution

When you want fewer voxels (e.g. to speed up analysis or plotting) but
don’t need arbitrarily defined target grids, use
[`downsample()`](https://bbuchsbaum.github.io/neuroim2/reference/downsample-methods.md)
on a `NeuroVec` or `NeuroVol`.

### 2.1 Downsample a 4D NeuroVec

``` r
# Downsample by a factor of 0.5 in each spatial dimension
vec_down_factor <- downsample(vec4, factor = 0.5)

dim(vec4)
#> [1] 64 64 25  4
dim(vec_down_factor)
#> [1] 32 32 12  4
spacing(vec4)
#> [1] 3.5 3.5 3.7
spacing(vec_down_factor)
#> [1] 7.00000 7.00000 7.70833
```

You can also specify a target spacing or output dimensions:

``` r
# Target spacing (mm)
vec_down_spacing <- downsample(vec4, spacing = c(4, 4, 4))

# Target spatial dimensions
vec_down_outdim  <- downsample(vec4, outdim = c(32, 32, 16))
```

Exactly one of `factor`, `spacing`, or `outdim` must be supplied. The
current implementation uses a simple box‑averaging scheme in space.

### 2.2 Downsample a 3D NeuroVol

The same interface applies to volumes:

``` r
vol_down_factor <- downsample(vol, factor = 0.5)

dim(vol)
#> [1] 64 64 25
dim(vol_down_factor)
#> [1] 32 32 12
```

This is useful for coarse previews or multi‑scale workflows.

## 3. Reorienting images

[`reorient()`](https://bbuchsbaum.github.io/neuroim2/reference/reorient-methods.md)
updates the mapping between voxel indices and physical coordinates
without changing the raw data array. This is helpful when you want a
canonical orientation (e.g. RAS) or need to match another dataset’s axis
directions.

For `NeuroSpace`, you supply three axis codes:

``` r
sp_lpi <- sp                      # assume input is LPI‑like
sp_ras <- reorient(sp_lpi, c("R", "A", "S"))

sp_lpi
#> <NeuroSpace> [3D] 
#> ── Geometry ──────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 64 x 64 x 25
#>   Spacing       : 3.5 x 3.5 x 3.7 mm
#>   Origin        : 112, -108, -46.2
#>   Orientation   : LAS
#>   Voxels        : 102,400
sp_ras
#> <NeuroSpace> [3D] 
#> ── Geometry ──────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 64 x 64 x 25
#>   Spacing       : 3.5 x 3.5 x 3.7 mm
#>   Origin        : -112, 108, 46.2
#>   Orientation   : RPI
#>   Voxels        : 102,400
```

For `NeuroVol` / `NeuroVec`, you typically reorient via their space:

``` r
vol_ras <- NeuroVol(as.array(vol), sp_ras)
dim(vol_ras)
#> [1] 64 64 25
spacing(vol_ras)
#> [1] 3.5 3.5 3.7
```

Note that
[`reorient()`](https://bbuchsbaum.github.io/neuroim2/reference/reorient-methods.md)
preserves physical positions; it only changes how grid indices are
interpreted in space. Coordinate transforms (`coord_to_grid`,
`grid_to_coord`) automatically respect the new orientation.

## 4. Putting it together

Common combinations:

- Use
  [`resample_to()`](https://bbuchsbaum.github.io/neuroim2/reference/resample_to.md)
  when you need to match a specific template grid (e.g. MNI space or
  another subject’s image).
- Use
  [`downsample()`](https://bbuchsbaum.github.io/neuroim2/reference/downsample-methods.md)
  when you just want fewer voxels while preserving overall FOV and
  orientation.
- Use
  [`reorient()`](https://bbuchsbaum.github.io/neuroim2/reference/reorient-methods.md)
  when you want a consistent anatomical convention (RAS/LPI) across
  images without resampling.
