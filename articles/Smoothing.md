# Smoothing and filtering with neuroim2

This vignette gives a practical overview of the main spatial and
spatio‑temporal smoothing tools in neuroim2:

- 3D spatial filters on `NeuroVol`:
  - [`gaussian_blur()`](https://bbuchsbaum.github.io/neuroim2/reference/gaussian_blur.md)
    — isotropic Gaussian smoothing
  - [`guided_filter()`](https://bbuchsbaum.github.io/neuroim2/reference/guided_filter.md)
    — edge‑preserving smoothing
  - [`bilateral_filter()`](https://bbuchsbaum.github.io/neuroim2/reference/bilateral_filter.md)
    — intensity‑aware spatial smoothing
  - [`laplace_enhance()`](https://bbuchsbaum.github.io/neuroim2/reference/laplace_enhance.md)
    — multi‑scale detail enhancement
- 4D filters on `NeuroVec`:
  - [`bilateral_filter_4d()`](https://bbuchsbaum.github.io/neuroim2/reference/bilateral_filter_4d.md)
    — joint spatial + temporal bilateral filter
  - [`cgb_filter()`](https://bbuchsbaum.github.io/neuroim2/reference/cgb_filter.md)
    — correlation‑guided bilateral smoothing based on a graph

We’ll work with the demo volume used elsewhere in the package. All code
assumes this setup:

``` r
set.seed(1)

demo_path <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
vol3d <- read_vol(demo_path)          # 3D volume
vec4d <- read_vec(demo_path)          # 4D NeuroVec (same file)

sp3  <- space(vol3d)
dims <- dim(vol3d)
```

## 1. Gaussian smoothing with `gaussian_blur()`

[`gaussian_blur()`](https://bbuchsbaum.github.io/neuroim2/reference/gaussian_blur.md)
applies an isotropic Gaussian kernel within an optional mask. It is
useful as a simple, fast baseline smoother.

- `sigma` (mm) controls how quickly weights decay with distance.
- `window` (voxels) sets the discrete kernel support (window = 1 →
  3×3×3).

``` r
blur_light <- gaussian_blur(vol3d, vol3d, sigma = 2, window = 1)
blur_strong <- gaussian_blur(vol3d, vol3d, sigma = 4, window = 2)

dim(blur_light)
#> [1] 64 64 25
```

Use smaller `sigma` and `window` for gentle smoothing; larger values
blur more but can wash out small structures.

## 2. Edge‑preserving guided filter with `guided_filter()`

[`guided_filter()`](https://bbuchsbaum.github.io/neuroim2/reference/guided_filter.md)
smooths while preserving edges by fitting local linear models between
the input and output.

- `radius` controls neighborhood size (in voxels).
- `epsilon` controls how strongly the filter smooths vs. preserves
  contrast.

``` r
gf_vol <- guided_filter(vol3d, radius = 4, epsilon = 0.7^2)
gf_vol
#> <DenseNeuroVol> [806.6 Kb] 
#> ── Spatial ───────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 64 x 64 x 25
#>   Spacing       : 3.5 x 3.5 x 3.7 mm
#>   Origin        : 112, -108, -46.2
#>   Orientation   : LAS
#> ── Data ──────────────────────────────────────────────────────────────────────── 
#>   Range         : [0.000, 1.000]
```

Compared to Gaussian smoothing, guided filtering tends to retain sharp
boundaries between tissues while denoising within regions.

## 3. Bilateral spatial filter with `bilateral_filter()`

[`bilateral_filter()`](https://bbuchsbaum.github.io/neuroim2/reference/bilateral_filter.md)
combines spatial distance and intensity similarity, reducing noise while
respecting edges.

``` r
bf_vol <- bilateral_filter(
  vol3d,
  spatial_sigma   = 2,
  intensity_sigma = 1,
  window          = 1
)
bf_vol
#> <DenseNeuroVol> [806.6 Kb] 
#> ── Spatial ───────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 64 x 64 x 25
#>   Spacing       : 3.5 x 3.5 x 3.7 mm
#>   Origin        : 112, -108, -46.2
#>   Orientation   : LAS
#> ── Data ──────────────────────────────────────────────────────────────────────── 
#>   Range         : [0.000, 0.990]
```

Key parameters:

- `spatial_sigma`: spatial scale (in mm).
- `intensity_sigma`: intensity scale (relative to global σ_I).
- `window`: discrete support (see
  [`?bilateral_filter`](https://bbuchsbaum.github.io/neuroim2/reference/bilateral_filter.md)).

Smaller `intensity_sigma` better preserves edges; larger values behave
more like pure Gaussian smoothing.

## 4. Laplacian enhancement with `laplace_enhance()`

[`laplace_enhance()`](https://bbuchsbaum.github.io/neuroim2/reference/laplace_enhance.md)
is designed for *sharpening* rather than smoothing. It uses a
multi‑layer, non‑local Laplacian scheme to enhance details.

``` r
sharp_vol <- laplace_enhance(vol3d, k = 2, patch_size = 3,
                             search_radius = 2, h = 0.7)
sharp_vol
#> <DenseNeuroVol> [806.6 Kb] 
#> ── Spatial ───────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 64 x 64 x 25
#>   Spacing       : 3.5 x 3.5 x 3.7 mm
#>   Origin        : 112, -108, -46.2
#>   Orientation   : LAS
#> ── Data ──────────────────────────────────────────────────────────────────────── 
#>   Range         : [-0.017, 1.013]
```

Use higher `h` or more layers `k` for stronger enhancement, but be
cautious about amplifying noise.

## 5. 4D bilateral smoothing with `bilateral_filter_4d()`

For 4D time‑series data (`NeuroVec`),
[`bilateral_filter_4d()`](https://bbuchsbaum.github.io/neuroim2/reference/bilateral_filter_4d.md)
extends the bilateral idea across space and time:

- `spatial_sigma`, `spatial_window` — as above.
- `intensity_sigma` — intensity similarity.
- `temporal_sigma`, `temporal_window` — temporal smoothing scale.
- `temporal_spacing` — units of the time axis (e.g., TR in seconds).

``` r
mask3d <- read_vol(system.file("extdata", "global_mask_v4.nii",
                               package = "neuroim2"))

bf4d <- bilateral_filter_4d(
  vec4d, mask3d,
  spatial_sigma   = 2,
  intensity_sigma = 1,
  temporal_sigma  = 1,
  spatial_window  = 1,
  temporal_window = 1,
  temporal_spacing = 1
)

dim(bf4d)
#> [1] 64 64 25  4
```

This is useful when you want joint spatial + temporal denoising that
still respects intensity boundaries (e.g. fMRI or 4D structural
sequences).

## 6. Correlation‑guided bilateral filtering with `cgb_filter()`

[`cgb_filter()`](https://bbuchsbaum.github.io/neuroim2/reference/cgb_filter.md)
implements correlation‑guided bilateral (CGB) smoothing: it builds a
sparse graph based on spatial distance and time‑series correlations,
then diffuses data over that graph.

Basic usage mirrors the bilateral interface:

``` r
cgf <- cgb_filter(
  vec4d,
  mask          = mask3d,
  spatial_sigma = 3,
  window        = NULL,      # auto from spatial_sigma and spacing
  corr_map      = "power",
  corr_param    = 2,
  topk          = 16,
  passes        = 1,
  lambda        = 1
)

dim(cgf)
#> [1] 64 64 25  4
```

Parameter intuition:

- `spatial_sigma` / `window`: spatial kernel scale and support (similar
  to Gaussian/bilateral).
- `corr_map`, `corr_param`: map pooled correlations to affinities; see
  [`?cgb_make_graph`](https://bbuchsbaum.github.io/neuroim2/reference/cgb_make_graph.md)
  for details.
- `topk`: max neighbors per voxel (sparsity vs. mixing).
- `passes`, `lambda`: diffusion strength (multiple passes or λ \< 1
  increase/temper smoothing).

Internally this calls
[`cgb_make_graph()`](https://bbuchsbaum.github.io/neuroim2/reference/cgb_make_graph.md)
to build the graph and then
[`cgb_smooth()`](https://bbuchsbaum.github.io/neuroim2/reference/cgb_smooth.md)
to diffuse over it. Use `return_graph = TRUE` to keep the graph for
reuse.

``` r
cg_out <- cgb_filter(vec4d, mask3d,
                     spatial_sigma = 3, window = NULL,
                     topk = 16, return_graph = TRUE)

str(cg_out$graph)
#> List of 5
#>  $ row_ptr : int [1:29533] 0 1 2 3 4 5 6 7 8 9 ...
#>  $ col_ind : int [1:29532] 474 475 476 477 480 481 482 483 536 537 ...
#>  $ val     : num [1:29532] 1 1 1 1 1 1 1 1 1 1 ...
#>  $ dims3d  : int [1:3] 64 64 25
#>  $ mask_idx: int [1:29532] 475 476 477 478 481 482 483 484 537 538 ...
```

## Choosing a smoother

- Use
  [`gaussian_blur()`](https://bbuchsbaum.github.io/neuroim2/reference/gaussian_blur.md)
  for simple, fast baseline smoothing.
- Use
  [`bilateral_filter()`](https://bbuchsbaum.github.io/neuroim2/reference/bilateral_filter.md)
  /
  [`bilateral_filter_4d()`](https://bbuchsbaum.github.io/neuroim2/reference/bilateral_filter_4d.md)
  when you want intensity‑aware, edge‑respecting smoothing.
- Use
  [`guided_filter()`](https://bbuchsbaum.github.io/neuroim2/reference/guided_filter.md)
  when you need edge preservation but want a simpler parameterization
  than full bilateral filters.
- Use
  [`laplace_enhance()`](https://bbuchsbaum.github.io/neuroim2/reference/laplace_enhance.md)
  when you want to *sharpen* structures rather than blur them.
- Use
  [`cgb_filter()`](https://bbuchsbaum.github.io/neuroim2/reference/cgb_filter.md)
  when similarity should be driven by *correlation of time‑series*
  (e.g., fMRI connectivity‑style smoothing) rather than raw intensity
  alone.
