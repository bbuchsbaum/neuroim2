# Extract Element from NeuroVecSeq

Extracts a single volume from a NeuroVecSeq object at the specified time
point.

## Usage

``` r
# S4 method for class 'NeuroVecSeq,numeric'
x[[i]]
```

## Arguments

- x:

  A NeuroVecSeq object

- i:

  Numeric index specifying the time point to extract

## Value

A NeuroVol object representing the extracted volume
