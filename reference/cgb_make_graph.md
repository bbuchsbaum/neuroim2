# Build a correlation-guided bilateral (CGB) graph

Computes a sparse row-stochastic graph whose weights combine spatial
proximity and pooled local time-series correlations. Supports optional
censoring weights, nuisance regression via weighted QR projectors,
leave-one-run-out graph construction, and robust down-weighting of
high-DVARS volumes.

## Usage

``` r
cgb_make_graph(
  runs,
  mask = NULL,
  window = 1L,
  spatial_sigma = 2,
  corr_map = c("power", "exp", "soft"),
  corr_param = 2,
  topk = 16L,
  leave_one_out = FALSE,
  run_weights = NULL,
  add_self = TRUE,
  time_weights = NULL,
  confounds = NULL,
  robust = c("none", "huber", "tukey"),
  robust_c = 1.345
)
```

## Arguments

- runs:

  A
  [`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md)
  or a list of `NeuroVec` objects (typically one per run).

- mask:

  Optional
  [`LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.md)/`NeuroVol`
  or logical array defining in-mask voxels. Defaults to all in-mask
  voxels.

- window:

  Integer half-width of the cubic spatial neighborhood (e.g., `1` yields
  a 3x3x3 window).

- spatial_sigma:

  Spatial Gaussian sigma in mm.

- corr_map:

  Mapping from pooled correlation to affinity; one of `"power"`,
  `"exp"`, or `"soft"`. The `"power"` and `"soft"` mappings rectify
  negative correlations, whereas `"exp"` preserves them (useful for
  sharpening more than smoothing).

- corr_param:

  Parameter for the chosen `corr_map` (gamma, tau, or r0 respectively).

- topk:

  Keep the strongest `k` neighbors after masking (0 keeps all).

- leave_one_out:

  Logical; if `TRUE` and multiple runs are provided, returns a list of
  graphs where run `u` excludes its own correlations.

- run_weights:

  Optional numeric weights per run used in Fisher-z pooling. Defaults to
  \\n_k - 3\\ (usable frames minus three) when omitted.

- add_self:

  Logical; always inject a tiny self-edge before normalization.

- time_weights:

  Optional list (or single vector) of per-run time weights \\w_t \in
  \[0,1\]\\ applied before correlation estimation. An intercept is
  always included so correlations are computed on weighted, demeaned
  series.

- confounds:

  Optional list (or single matrix) of per-run confound regressors to
  project out prior to correlation estimation.

- robust:

  One of `"none"`, `"huber"`, or `"tukey"`; when not `"none"` an
  additional DVARS-style reweighting is applied.

- robust_c:

  Tuning constant for the robust weights (Huber/Tukey).

## Value

A list containing `row_ptr`, `col_ind`, `val`, `dims3d`, and `mask_idx`,
or (if `leave_one_out=TRUE`) a list of such graphs.

## Details

Graph construction overview: - Neighborhood: For each in-mask voxel i,
consider a cubic spatial window of half-width `window` (i.e.,
(2\*window+1)^3 candidates). Candidates outside the mask or bounds are
ignored. - Spatial kernel: For a candidate j at physical distance d_ij
(mm), assign a spatial weight w_s = exp(-d_ij^2 / (2 \*
spatial_sigma^2)). Distances use `spacing(spatial_space)`. - Correlation
pooling: Compute Pearson correlation r_k(i,j) within each run k
(optionally after nuisance projection/weights), transform to Fisher-z,
pool across runs with weights \\\omega_k\\ (default \\n_k - 3\\), then
back-transform to r_pool via tanh. - Correlation-to-affinity mapping
(`corr_map`): \* `"power"` (mode=0): a(r) = r^gamma for r\>0 else 0.
Parameter = gamma. \* `"exp"` (mode=1): a(r) = exp(- (1 - r)^2 / (2 \*
tau^2)). Parameter = tau. \* `"soft"` (mode=2): a(r) = max(r - r0, 0).
Parameter = r0. - Combined weight: w_ij = w_s(i,j) \* a(r_pool(i,j)). If
`topk > 0`, keep the strongest `topk` neighbors. Optionally inject a
small self-edge when `add_self=TRUE`. Finally, row-normalize to obtain a
stochastic W.

Parameter guidance: - `spatial_sigma` is in mm. A typical choice is 1-2x
the voxel size (e.g., 2-4 mm for 2 mm isotropic). Larger values increase
spatial mixing. - `window` controls support; a good rule is
`ceiling(2 * spatial_sigma / min(spacing))`. - `corr_map`: use `"power"`
with `corr_param = 2` for robust smoothing; `"exp"` with `tau ~ 0.5-1.5`
retains sign information; `"soft"` with `r0 ~ 0.1-0.3` thresholds weak
correlations. - `topk`: 8-32 is a practical range; higher values densify
the graph and increase compute/memory. - `leave_one_out`: for multi-run
inputs, enabling this prevents a run from using its own correlations
when building its graph (mitigates leakage).

Nuisance/time weights/robust options: - If `confounds`/`time_weights`
specified (or `robust != "none"`), per-run weighted QR projectors are
used; correlations are computed on the projected, weighted series.
Robust options add DVARS-like down-weighting. - If the only request is
an implicit intercept (no actual confounds, no time weights, no robust),
the baseline builder is used for identical results.

Output and usage: - Returns CSR arrays (`row_ptr`, `col_ind`, `val`)
that define a row-stochastic matrix W over masked voxels. Apply with
`cgb_smooth`. - Complexity scales with number of masked voxels times
neighborhood size (limited by `topk`). Memory proportional to number of
retained edges.

## Examples

``` r
# \donttest{
vec <- read_vec(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))
mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))
# Build a graph with spatial sigma in mm and keep 16 neighbors
G <- cgb_make_graph(vec, mask, spatial_sigma = 3, window = 2, topk = 16)
# Smooth with one pass (pure diffusion)
sm <- cgb_smooth(vec, G, passes = 1, lambda = 1)
# }
```
