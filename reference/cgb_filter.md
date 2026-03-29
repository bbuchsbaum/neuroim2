# Correlation-guided bilateral filtering (convenience wrapper)

High-level interface that builds a correlation-guided bilateral (CGB)
graph with sensible defaults (similar to the bilateral filter interface)
and immediately applies it to smooth the data.

## Usage

``` r
cgb_filter(
  runs,
  mask = NULL,
  spatial_sigma = 2,
  window = NULL,
  corr_map = c("power", "exp", "soft"),
  corr_param = 2,
  topk = 16L,
  passes = 1L,
  lambda = 1,
  leave_one_out = FALSE,
  run_weights = NULL,
  add_self = TRUE,
  time_weights = NULL,
  confounds = NULL,
  robust = c("none", "huber", "tukey"),
  robust_c = 1.345,
  return_graph = FALSE
)
```

## Arguments

- runs:

  A
  [`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md)
  or a list of `NeuroVec`.

- mask:

  Optional
  [`LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.md)/`NeuroVol`
  or logical array for spatial masking. Defaults to in-mask voxels.

- spatial_sigma:

  Spatial Gaussian sigma in mm. Used both for weighting and, when
  `window` is `NULL`, to auto-choose the neighborhood size.

- window:

  Integer half-width of the cubic neighborhood. If `NULL`, it is
  computed as `ceiling(2 * spatial_sigma / min(spacing))` and at least
  1.

- corr_map:

  Mapping from pooled correlation to affinity; one of `"power"`,
  `"exp"`, or `"soft"`. Defaults to `"power"`.

- corr_param:

  Parameter for `corr_map` (gamma/tau/r0 respectively).

- topk:

  Keep strongest `k` neighbors (0 keeps all). Defaults to 16.

- passes:

  Number of smoothing passes (\>=1). Defaults to 1.

- lambda:

  Blend factor in \[0,1\] per pass. Defaults to 1 (pure diffusion).

- leave_one_out:

  If `TRUE` and multiple runs are supplied, builds LORO graphs and
  returns a list of smoothed runs.

- run_weights:

  Optional numeric weights per run for Fisher-z pooling.

- add_self:

  Logical; add a tiny self-edge before normalization.

- time_weights:

  Optional list (or single vector) of per-run time weights.

- confounds:

  Optional list (or single matrix) of per-run confounds.

- robust:

  One of `"none"`, `"huber"`, or `"tukey"`.

- robust_c:

  Tuning constant for robust weights.

- return_graph:

  Logical; if `TRUE`, also return the graph(s) alongside the smoothed
  data.

## Value

If `leave_one_out=FALSE`, a smoothed `NeuroVec`. If
`leave_one_out=TRUE`, a list of smoothed `NeuroVec`. When
`return_graph=TRUE`, returns a list with elements `result` and `graph`
(single object or lists accordingly).

## Details

This is a convenience front-end to `cgb_make_graph` and `cgb_smooth`
with a bilateral-like interface: - If `window` is `NULL`, it is chosen
as `ceiling(2 * spatial_sigma / min(spacing))` (at least 1). Larger
windows allow correlations over more distant neighbors, at the cost of
extra compute and memory. - `spatial_sigma` (in mm) controls how quickly
spatial weights fall with distance. Small values emphasize very local
structure; larger values mix information over a wider spatial
footprint. - `corr_map` and `corr_param` set how pooled correlations are
turned into edge weights: \* `"power"`: `a(r) = r^gamma` for `r > 0`.
Larger `gamma` (e.g., 3-4) strongly emphasizes high correlations and
produces more edge-preserving, patchy smoothing; smaller values (e.g.,
1-2) behave more like standard correlation-weighted smoothing. \*
`"exp"`: Gaussian on `1 - r` with scale `tau`. Small `tau` keeps only
very similar time-series; larger `tau` makes the filter closer to a
spatial Gaussian while still respecting sign. \* `"soft"`:
`a(r) = max(r - r0, 0)`. Increasing `r0` discards more weak correlations
and tends to sharpen edges but can make the result more
piecewise-constant. - `topk` limits each voxel to at most `k` strongest
neighbors. Smaller `topk` yields sparser, more anisotropic graphs
(cheaper but sometimes less smooth); larger `topk` increases mixing and
memory. - `passes` and `lambda` control diffusion strength. With
`lambda = 1`, each pass applies pure graph diffusion; multiple passes
compound smoothing. Choosing `lambda < 1` blends each pass with the
identity and can prevent over-smoothing when using more passes. -
Setting `leave_one_out = TRUE` for multi-run inputs builds a separate
graph for each run that excludes its own correlations, which reduces
information leakage in cross-validation or decoding workflows. -
`time_weights`, `confounds`, and `robust`/`robust_c` adjust how
time-points contribute to the correlation estimates. Down- weighting
high-motion/high-DVARS frames (via `make_time_weights` and
`robust != "none"`) will typically yield smoother, less noisy graphs but
can also reduce effective temporal degrees of freedom. - Use
`return_graph = TRUE` when you plan to reuse the constructed graph(s)
with `cgb_smooth` or inspect their sparsity pattern.

## Examples

``` r
# \donttest{
vec <- read_vec(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))
mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))

# Auto window from spatial_sigma and spacing, single pass
out <- cgb_filter(vec, mask, spatial_sigma = 3, window = NULL, topk = 16)

# Stronger diffusion with two passes and lambda < 1
out2 <- cgb_filter(vec, mask, spatial_sigma = 4, window = NULL,
                   passes = 2, lambda = 0.7)
# }
```
