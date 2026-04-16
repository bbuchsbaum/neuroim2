# Generic function to extract a sub-vector from a `NeuroVec` object.

Extracts a subset of volumes from a file-backed neuroimaging vector and
returns them as a dense (in-memory) vector.

Extracts a subsequence of volumes from a NeuroVecSeq object.

## Usage

``` r
sub_vector(x, i, ...)

# S4 method for class 'FileBackedNeuroVec,numeric'
sub_vector(x, i)

# S4 method for class 'NeuroVec,numeric'
sub_vector(x, i)

# S4 method for class 'NeuroVec,character'
sub_vector(x, i)

# S4 method for class 'NeuroVecSeq,numeric'
sub_vector(x, i)

# S4 method for class 'NeuroVecSeq,numeric'
sub_vector(x, i)

# S4 method for class 'AbstractSparseNeuroVec,numeric'
sub_vector(x, i)
```

## Arguments

- x:

  A NeuroVecSeq object

- i:

  Numeric vector of indices specifying the time points to extract

- ...:

  additional arguments

## Value

A `NeuroVec` object that is a sub-sequence of the supplied object.

A NeuroVecSeq object containing the extracted subsequence

## Details

This method efficiently reads only the requested volumes from disk,
converting them to an in-memory representation. The spatial metadata is
preserved but adjusted to reflect the new number of volumes.

Memory usage is proportional to the number of volumes requested, not the
size of the full dataset.

## Examples

``` r
bvec <- NeuroVec(array(rnorm(24*24*24*24), c(24,24,24,24)), NeuroSpace(c(24,24,24,24), c(1,1,1)))
vec <- sub_vector(bvec,1:2)
all.equal(2, dim(vec)[4])
#> [1] TRUE

vec <- sub_vector(bvec, c(1,3,5,7))
all.equal(4, dim(vec)[4])
#> [1] TRUE

mask <- LogicalNeuroVol(rep(TRUE, 24*24*24), NeuroSpace(c(24,24,24), c(1,1,1)))
svec <- SparseNeuroVec(array(rnorm(24*24*24*24), c(24,24,24,24)),
NeuroSpace(c(24,24,24,24), c(1,1,1)), mask)
vec <- sub_vector(svec, c(1,3,5))
all.equal(3, dim(vec)[4])
#> [1] TRUE
```
