# Cut a vector-valued object into a list of sub-blocks

Splits a vector-valued object into a list of sub-blocks defined by a
vector of indices.

## Usage

``` r
split_blocks(x, indices, ...)

# S4 method for class 'NeuroVec,integer'
split_blocks(x, indices, ...)

# S4 method for class 'NeuroVec,factor'
split_blocks(x, indices, ...)

# S4 method for class 'NeuroVec,factor'
split_blocks(x, indices, ...)
```

## Arguments

- x:

  a vector-valued object

- indices:

  a vector of indices defining the sub-blocks. Must match the length of
  the input vector.

- ...:

  additional arguments

## Value

A `list` of sub-blocks, where each sub-block contains the elements from
`x` corresponding to the matching `indices`.

## Examples

``` r
# Create a 4D neuroimaging vector with 20 timepoints
space <- NeuroSpace(c(10,10,10,20), c(1,1,1))
vec <- NeuroVec(array(rnorm(10*10*10*20), c(10,10,10,20)), space)

# Split into 4 blocks by assigning timepoints to blocks 1-4 repeatedly
block_indices <- rep(1:4, length.out=20)
blocks <- split_blocks(vec, block_indices)
```
