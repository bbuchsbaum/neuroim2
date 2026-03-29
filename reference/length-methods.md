# Get length of NeuroVec object

Returns the number of time points (4th dimension) in a NeuroVec object.
This represents the temporal dimension of the neuroimaging data.

Returns the total number of time points across all vectors in the
sequence

## Usage

``` r
# S4 method for class 'ClusteredNeuroVec'
length(x)

# S4 method for class 'NeuroVec'
length(x)

# S4 method for class 'NeuroVecSeq'
length(x)

# S4 method for class 'ROIVol'
length(x)

# S4 method for class 'ROICoords'
length(x)
```

## Arguments

- x:

  A NeuroVecSeq object

## Value

Integer length (total number of time points)

An integer representing the number of coordinates in the ROICoords
object.
