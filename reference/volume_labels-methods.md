# Get per-volume labels for a `NeuroVec`.

Get per-volume labels for a `NeuroVec`.

## Usage

``` r
volume_labels(x)

# S4 method for class 'NeuroVec'
volume_labels(x)
```

## Arguments

- x:

  A `NeuroVec` or compatible object.

## Value

A character vector of labels. Returns `character(0)` when no labels are
defined.

A character vector of per-volume labels, or `character(0)` if none are
set.

## Examples

``` r
sp <- NeuroSpace(c(2, 2, 2, 3), c(1, 1, 1))
vec <- NeuroVec(array(1:24, dim = c(2, 2, 2, 3)), sp,
                volume_labels = c("baseline", "task", "rest"))
volume_labels(vec)
#> [1] "baseline" "task"     "rest"    
```
