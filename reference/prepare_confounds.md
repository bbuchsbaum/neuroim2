# Prepare weighted nuisance projectors for each run

Converts per-run confound matrices and time weights into orthonormal
projectors that can be consumed directly by the C++ graph builder. Each
run produces \\Q_k\\ (columns spanning the weighted confound space) and
\\\sqrt{w_k}\\ (per-time-point square-root weights). Supplying a `NULL`
confound matrix yields a zero-column projector, enabling pure
time-weighting without regression.

## Usage

``` r
prepare_confounds(
  confounds,
  time_weights = NULL,
  run_lengths,
  include_intercept = TRUE,
  center_cols = TRUE,
  scale_cols = FALSE
)
```

## Arguments

- confounds:

  List of matrices (\\T_k \times p_k\\), or a single matrix reused for
  all runs. Each row corresponds to a time point.

- time_weights:

  Optional list of numeric vectors (length \\T_k\\) or a single vector
  reused for every run. If `NULL`, unit weights are used.

- run_lengths:

  Integer vector with the number of time points per run. Required when
  any run has both `confounds=NULL` and `time_weights=NULL`.

- include_intercept:

  Logical; prepend a column of ones before QR (default TRUE).

- center_cols:

  Logical; center each confound column before weighting (default TRUE).

- scale_cols:

  Logical; scale columns to unit variance before weighting (default
  FALSE).

## Value

A list with elements `Q_list` (list of matrices) and `sqrtw_list` (list
of numeric vectors). Each entry has the same length as `run_lengths`.
