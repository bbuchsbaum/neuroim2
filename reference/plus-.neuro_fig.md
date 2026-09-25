# Arithmetic on device-fitted figures

`+` and `&` work as for any patchwork; the additions are replayed when
the figure is re-fitted to the device it is drawn on.

## Usage

``` r
# S3 method for class 'neuro_fig'
e1 + e2

# S3 method for class 'neuro_fig'
e1 & e2

# S3 method for class 'neuro_fig'
e1 | e2

# S3 method for class 'neuro_fig'
e1/e2
```

## Arguments

- e1:

  A figure returned by a `plot_*` function.

- e2:

  A ggplot2 or patchwork component.

## Value

A figure of the same class.
