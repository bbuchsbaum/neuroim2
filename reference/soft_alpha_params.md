# Self-tuning nonlinear opacity curve for statistical overlays

Parameters of the curve used by `plot_overlay(ov_alpha_mode = "soft")`.
For magnitude \\m\\, opacity is \$\$alpha(m) = f + (1 -
f)\\\mathrm{clamp}((m - lo) / (hi - lo), 0, 1)^{\gamma},\$\$ where \\f\\
is `alpha_floor`. Plotting then hides values below the hard threshold
and multiplies by `ov_alpha`.

## Usage

``` r
soft_alpha_params(
  mags,
  thresh = 0,
  cap = NULL,
  gamma = NULL,
  alpha_mid = 0.2,
  gamma_min = 1.5,
  gamma_max = 5,
  knee = NULL,
  alpha_floor = 0
)
```

## Arguments

- mags:

  Numeric vector of overlay magnitudes (typically `abs(values)`).

- thresh:

  Hard threshold; used as the knee when `> 0` and `knee` is not given.

- cap:

  Optional upper magnitude anchor (opacity 1).

- gamma:

  Optional fixed exponent; `NULL` auto-tunes it.

- alpha_mid:

  Target opacity for the median supra-knee magnitude when `gamma` is
  tuned.

- gamma_min, gamma_max:

  Clamp range for the tuned exponent.

- knee:

  Optional non-negative lower magnitude anchor, overriding the
  threshold/median policy. Use `0` to ramp from zero.

- alpha_floor:

  Minimum opacity (0 to 1) above the knee, before `ov_alpha` and the
  hard threshold are applied.

## Value

A list with `lo`, `hi`, `gamma` and `alpha_floor`.

## Details

The knee `lo` defaults to the threshold or, when no threshold is set, to
the median non-zero magnitude (a robust noise-floor proxy). The cap `hi`
defaults to the largest magnitude. When `gamma` is not given it is tuned
so the median supra-knee magnitude maps to `alpha_mid` (before the
floor), and clamped to `[gamma_min, gamma_max]` so the default curve
stays convex; an explicit `gamma` bypasses the clamp. Fix `knee`, `cap`
and `gamma` to reuse one curve across datasets.

## See also

[`plot_overlay`](https://bbuchsbaum.github.io/neuroim2/reference/plot_overlay.md)

## Examples

``` r
p <- soft_alpha_params(0:8, knee = 0, cap = 3, gamma = 0.7, alpha_floor = 0.15)
m <- c(0, 1.6, 2.1, 3, 8)
p$alpha_floor + (1 - p$alpha_floor) *
  pmin(pmax((m - p$lo) / (p$hi - p$lo), 0), 1)^p$gamma
#> [1] 0.1500000 0.6974167 0.8121975 1.0000000 1.0000000
```
