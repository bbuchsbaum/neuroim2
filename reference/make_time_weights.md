# Build smooth time weights from motion/outlier metrics

Creates per-time-point weights \\w_t \in \[0, 1\]\\ by smoothly
combining framewise displacement (FD), DVARS, and spike/outlier scores.
Each series is transformed through a soft logistic ramp so that values
beyond the specified thresholds receive progressively lower weights
instead of hard 0/1 decisions.

## Usage

``` r
make_time_weights(
  fd = NULL,
  dvars = NULL,
  spike = NULL,
  fd_thr = 0.5,
  dvars_z = 2.5,
  spike_z = 5,
  fd_soft = 0.1,
  dvars_soft = 0.25,
  combine = c("min", "prod")
)
```

## Arguments

- fd:

  Optional numeric vector of framewise displacement values.

- dvars:

  Optional numeric vector of DVARS values.

- spike:

  Optional numeric vector with spike/outlier magnitudes.

- fd_thr:

  Threshold (in mm) where FD weights start to drop (default 0.5).

- dvars_z:

  Z-threshold applied to the standardized DVARS series (default 2.5).

- spike_z:

  Z-threshold applied to the standardized spike series (default 5).

- fd_soft:

  Logistic softness (in mm) controlling the slope around `fd_thr`.

- dvars_soft:

  Logistic softness for the DVARS z-scores.

- combine:

  Either `"min"` (take the minimum weight per TR) or `"prod"` (multiply
  all weights).

## Value

Numeric vector of weights in \\\[0,1\]\\ with length equal to the
provided series. At least one of `fd`, `dvars`, or `spike` must be
supplied.
