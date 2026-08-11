# Convert ClusteredNeuroVol to a base array

Ensures that clustered volumes dispatch through the \`as.array\` S4
generic and return dense arrays of cluster labels aligned to the
underlying space.

Provides an \`as.array\` S4 method so sparse volumes can be coerced with
the same syntax used for dense objects.

The S4 slots of an object that extends `array` are stored as attributes
on that array, so the default coercion handed back the voxel data still
carrying `space`, `header` and `class`. Anyone calling
[`as.array()`](https://bbuchsbaum.github.io/neuroim2/reference/as.array.md)
wants the numbers without the wrapper; these methods give exactly that.

## Usage

``` r
# S4 method for class 'ClusteredNeuroVol'
as.array(x, ...)

# S4 method for class 'DenseNeuroVec'
as.array(x, ...)

# S4 method for class 'SparseNeuroVol'
as.array(x, ...)

# S4 method for class 'DenseNeuroVol'
as.array(x, ...)

# S4 method for class 'NeuroSlice'
as.array(x, ...)

# S4 method for class 'SparseNeuroVec'
as.array(x, ...)
```

## Arguments

- x:

  A SparseNeuroVec object.

- ...:

  Additional arguments (currently ignored).

## Value

A dense array of cluster ids.

A dense array with voxel values at their spatial locations and zeros
elsewhere.

A base array with only a `dim` attribute.

A dense 4D array with sparse values inserted at mask indices and zeros
elsewhere.
