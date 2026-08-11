# Enhance an unsmoothed statistical map for visualization

De-speckles a noisy ("salt-and-pepper") statistical map so it renders
cleanly as an overlay, *without* the peak-depressing and
cluster-smearing that a plain Gaussian blur causes. This is a
display-oriented preprocessor: it is not intended to produce maps for
statistical inference (it deliberately sharpens for visual punch and
does not preserve the null distribution).

## Usage

``` r
enhance_stat_map(
  vol,
  mask,
  despike = TRUE,
  despike_k = 3.5,
  despike_window = 1L,
  method = c("guided", "bilateral", "gaussian"),
  radius = 2,
  epsilon = NULL,
  spatial_sigma = 2,
  intensity_sigma = 1,
  detail_gain = 1.5,
  verbose = FALSE
)
```

## Arguments

- vol:

  A
  [`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
  statistical map to enhance.

- mask:

  An optional
  [`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
  /
  [`LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.md)
  restricting the region processed. If omitted, all non-zero, finite
  voxels of `vol` are used.

- despike:

  Logical; run the selective median despike stage (default `TRUE`).

- despike_k:

  Robust-deviation threshold for flagging impulses. Lower is more
  aggressive (default 3.5).

- despike_window:

  Integer half-width of the despike neighbourhood (default 1L =\>
  3x3x3).

- method:

  Edge-preserving base filter: one of `"guided"` (default),
  `"bilateral"`, or `"gaussian"`.

- radius:

  Integer box radius for the base filter and signal-weight estimation
  (default 2).

- epsilon:

  Guided-filter regularisation. If `NULL` (default) it is set to the
  estimated noise variance. Only used when `method = "guided"`.

- spatial_sigma:

  Spatial Gaussian SD for `method = "bilateral"` or `"gaussian"`
  (default 2).

- intensity_sigma:

  Intensity Gaussian SD for `method = "bilateral"` (default 1).

- detail_gain:

  Strength of the signal-gated unsharp stage. `0` gives a pure denoise;
  `1` restores original detail in signal regions; values `> 1` sharpen
  for display (default 1.5). Set to 0 to disable.

- verbose:

  Logical; report the number of despiked voxels and the noise estimate
  (default `FALSE`).

## Value

A
[`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.md)
with the enhanced map (zero outside the mask).

## Details

The enhancement runs in three stages, each reusing the package's
existing filtering primitives:

1.  **Selective despike.** Isolated impulsive voxels (flagged by a
    robust deviation from their local mean) are replaced by the median
    of their in-mask, non-impulse neighbours. All other voxels keep
    their exact value, so true clusters are preserved.

2.  **Edge-preserving base smooth.** An edge-preserving filter
    (`"guided"` by default; see
    [`guided_filter`](https://bbuchsbaum.github.io/neuroim2/reference/guided_filter.md)
    /
    [`bilateral_filter`](https://bbuchsbaum.github.io/neuroim2/reference/bilateral_filter.md)
    /
    [`gaussian_blur`](https://bbuchsbaum.github.io/neuroim2/reference/gaussian_blur.md))
    removes residual low-amplitude grain. With the guided filter,
    in-cluster voxels (local variance \\\gg\\ noise variance) are
    returned at essentially full amplitude, so peaks are not depressed.

3.  **Signal-gated unsharp.** Local contrast is boosted by
    `detail_gain`, weighted per-voxel by a guided-style signal weight
    \\w = \mathrm{var}/(\mathrm{var} + \sigma^2\_{noise})\\. Because \\w
    \to 0\\ in flat/noisy regions, the despeckled background stays
    smooth while real clusters gain crisp, punchy edges.

The noise scale \\\sigma\_{noise}\\ is estimated robustly (MAD of the
high-frequency residual within the mask) and is used both for the
default guided `epsilon` and for the signal weight, so the method adapts
to the scale of the input map (e.g. a z-map versus a raw beta map).

## See also

[`guided_filter`](https://bbuchsbaum.github.io/neuroim2/reference/guided_filter.md),
[`bilateral_filter`](https://bbuchsbaum.github.io/neuroim2/reference/bilateral_filter.md),
[`gaussian_blur`](https://bbuchsbaum.github.io/neuroim2/reference/gaussian_blur.md),
[`plot_overlay`](https://bbuchsbaum.github.io/neuroim2/reference/plot_overlay.md),
[`plot_ortho`](https://bbuchsbaum.github.io/neuroim2/reference/plot_ortho.md)

## Examples

``` r
# \donttest{
mask <- read_vol(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))
# A noisy synthetic stat map on the mask grid
set.seed(1)
noisy <- mask * 0
idx <- which(mask > 0)
vals <- rnorm(length(idx))
vals[sample(length(idx), length(idx) %/% 20)] <- 8  # salt-and-pepper spikes
noisy[idx] <- vals
stat <- NeuroVol(noisy, space(mask))

clean <- enhance_stat_map(stat, mask)
# }
```
