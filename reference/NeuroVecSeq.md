# NeuroVecSeq: A Container for Sequential NeuroVec Objects

The NeuroVecSeq class provides a container for managing a sequence of
NeuroVec objects, particularly useful for handling time series or
multi-session neuroimaging data where each segment may have different
lengths.

Constructs a NeuroVecSeq object to represent a variable-length sequence
of NeuroVec objects. This is particularly useful for managing time
series data where different segments may have different lengths.

## Usage

``` r
NeuroVecSeq(...)
```

## Arguments

- ...:

  One or more instances of type
  [`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md).

## Value

A NeuroVecSeq object containing:

- The provided NeuroVec objects

- Associated space information

- Length information for each vector

## Details

NeuroVecSeq objects store:

- A list of NeuroVec objects, each potentially with different time
  dimensions

- The lengths of each constituent NeuroVec

- A combined NeuroSpace object representing the total space

The class provides methods for:

- Accessing individual time points across all vectors

- Extracting subsequences

- Computing statistics across the sequence

- Linear access to the underlying data

The function performs several validations:

- Ensures all inputs are NeuroVec objects

- Verifies spatial compatibility

- Combines spatial information appropriately

## Methods

- \[\[:

  Extract a single volume at a specified time point

- length:

  Get the total number of time points

- sub_vector:

  Extract a subsequence of volumes

- linear_access:

  Access data linearly across all vectors

## See also

[`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md)
for the base vector class,
[`NeuroSpace`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.md)
for spatial information

## Examples

``` r
# Create some example NeuroVec objects
v1 <- NeuroVec(array(0, c(5, 5, 5, 2)),
               space = NeuroSpace(dim = c(5, 5, 5, 2)))
v2 <- NeuroVec(array(1, c(5, 5, 5, 4)),
               space = NeuroSpace(dim = c(5, 5, 5, 4)))
v3 <- NeuroVec(array(2, c(5, 5, 5, 6)),
               space = NeuroSpace(dim = c(5, 5, 5, 6)))

# Combine them into a sequence
vs <- NeuroVecSeq(v1, v2, v3)

# Access properties
length(vs)  # Total time points
#> [1] 12
vs[[5]]     # Get the 5th volume
#> <DenseNeuroVol> [7.6 Kb] 
#> ── Spatial ───────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 5 x 5 x 5
#>   Spacing       : 1 x 1 x 1 mm
#>   Origin        : 0, 0, 0
#>   Orientation   : RAS
#> ── Data ──────────────────────────────────────────────────────────────────────── 
#>   Range         : [1.000, 1.000]

# Extract a subsequence
sub_seq <- sub_vector(vs, 1:5)


# Create sample vectors
v1 <- NeuroVec(array(0, c(5, 5, 5, 2)),
               space = NeuroSpace(dim = c(5, 5, 5, 2)))
v2 <- NeuroVec(array(0, c(5, 5, 5, 4)),
               space = NeuroSpace(dim = c(5, 5, 5, 4)))

# Combine into sequence
vs <- NeuroVecSeq(v1, v2)
print(vs)
#> <NeuroVecSeq> [2 vectors] 
#> ── Sequence ──────────────────────────────────────────────────────────────────── 
#>     [1]         : DenseNeuroVec 5x5x5 x 2t
#>     [2]         : DenseNeuroVec 5x5x5 x 4t

```
