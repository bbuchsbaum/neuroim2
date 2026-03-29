# as.logical

Convert NeuroVol to
[`LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.md)

## Usage

``` r
# S4 method for class 'NeuroVol'
as.logical(x)

# S4 method for class 'ROIVol'
as.logical(x)
```

## Arguments

- x:

  the object

## Value

an instance of
[`LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.md)

## Details

the image values will be converted to using R base function `as.logical`
and wrapped in `LogicalNeuroVol`
