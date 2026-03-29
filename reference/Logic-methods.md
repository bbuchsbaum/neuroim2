# Logic Operations for Neuroimaging Volumes

Methods for performing logical operations (`&` and `|`) on neuroimaging
volume objects. Results are always returned as
[`LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.md)
objects that preserve spatial metadata.

## Usage

``` r
# S4 method for class 'DenseNeuroVol,DenseNeuroVol'
Logic(e1, e2)

# S4 method for class 'SparseNeuroVol,SparseNeuroVol'
Logic(e1, e2)

# S4 method for class 'SparseNeuroVol,NeuroVol'
Logic(e1, e2)

# S4 method for class 'NeuroVol,SparseNeuroVol'
Logic(e1, e2)

# S4 method for class 'NeuroVol,logical'
Logic(e1, e2)

# S4 method for class 'logical,NeuroVol'
Logic(e1, e2)
```

## Arguments

- e1, e2:

  Neuroimaging volume objects or logical values.

## Value

A
[`LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.md).

## Examples

``` r
sp <- NeuroSpace(c(5L, 5L, 5L), c(1, 1, 1))
v1 <- DenseNeuroVol(array(sample(0:1, 125, replace = TRUE), c(5, 5, 5)), sp)
v2 <- DenseNeuroVol(array(sample(0:1, 125, replace = TRUE), c(5, 5, 5)), sp)
intersection <- v1 & v2
union_mask  <- v1 | v2
```
