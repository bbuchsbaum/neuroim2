# Convert SparseNeuroVol to a base vector

Supplies an \`as.vector\` S4 method that flattens sparse volumes to a
dense vector, keeping the same voxel ordering as \`as.array\`.

## Usage

``` r
# S4 method for class 'SparseNeuroVol'
as.vector(x, mode = "any")
```

## Arguments

- x:

  A \`SparseNeuroVol\` instance.

- mode:

  Optional coercion mode (see \[base::as.vector\]).

## Value

A vector of length \`prod(dim(x))\`.
