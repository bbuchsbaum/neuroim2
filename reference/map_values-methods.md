# Map Values from One Set to Another Using a User-supplied Lookup Table

This function maps values from one set to another using a lookup table
provided by the user.

## Usage

``` r
map_values(x, lookup)

# S4 method for class 'NeuroVol,list'
map_values(x, lookup)

# S4 method for class 'NeuroVol,matrix'
map_values(x, lookup)
```

## Arguments

- x:

  The object from which values will be mapped.

- lookup:

  The lookup table. The first column is the "key" and the second column
  is the "value".

## Value

An object of the same class as `x`, in which the original values have
been replaced with the lookup table values.

## Examples

``` r
x <- NeuroSpace(c(10, 10, 10), c(1, 1, 1))
vol <- NeuroVol(sample(1:10, 10 * 10 * 10, replace = TRUE), x)

## Lookup table is a list
lookup <- lapply(1:10, function(i) i * 10)
names(lookup) <- 1:10
ovol <- map_values(vol, lookup)

## Lookup table is a matrix. The first column is the key, and the second column is the value
names(lookup) <- 1:length(lookup)
lookup.mat <- cbind(as.numeric(names(lookup)), unlist(lookup))
ovol2 <- map_values(vol, lookup.mat)
all.equal(as.vector(ovol2), as.vector(ovol))
#> [1] TRUE
```
