# Fill Disjoint Sets of Values with the Output of a Function

This function splits an object into disjoint sets of values based on a
factor, applies a specified function to each set, and returns a new
object with the original values replaced by the function's output.

## Usage

``` r
split_fill(x, fac, FUN)

# S4 method for class 'NeuroVol,factor,function'
split_fill(x, fac, FUN)
```

## Arguments

- x:

  The object to split.

- fac:

  The `factor` to split by.

- FUN:

  The function used to summarize the sets.

## Value

An object of the same class as `x`, with values replaced by the output
of `FUN`.

## Details

The `FUN` function can either return a scalar for each input vector or a
vector equal to the length of the input vector. If it returns a scalar,
every voxel in the set will be filled with that value in the output
vector.

## Examples

``` r
## Summarize with mean -- FUN returns a scalar
x <- NeuroSpace(c(10, 10, 10), c(1, 1, 1))
vol <- NeuroVol(rnorm(10 * 10 * 10), x)
fac <- factor(rep(1:10, length.out=1000))
ovol.mean <- split_fill(vol, fac, mean)
identical(dim(ovol.mean), dim(vol))
#> [1] TRUE
length(unique(as.vector(ovol.mean))) == 10
#> [1] TRUE

## Transform by reversing vector -- FUN returns a vector
ovol2 <- split_fill(vol, fac, rev)
```
