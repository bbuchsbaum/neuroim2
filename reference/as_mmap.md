# Convert a NeuroVec to a memory-mapped representation

Generic for converting neuroimaging vectors to a memory-mapped
[`MappedNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/MappedNeuroVec-class.md)
on disk (when possible).

Methods for the `as_mmap` generic, which convert various neuroimaging
vector types to a
[`MappedNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/MappedNeuroVec-class.md)
backed by an on-disk NIfTI file.

## Usage

``` r
as_mmap(x, file = NULL, ...)

# S4 method for class 'MappedNeuroVec'
as_mmap(x, file = NULL, ...)

# S4 method for class 'FileBackedNeuroVec'
as_mmap(x, file = NULL, ...)

# S4 method for class 'NeuroVec'
as_mmap(x, file = NULL, data_type = "FLOAT", overwrite = FALSE, ...)

# S4 method for class 'SparseNeuroVec'
as_mmap(x, file = NULL, data_type = "FLOAT", overwrite = FALSE, ...)
```

## Arguments

- x:

  A neuroimaging vector (`NeuroVec`, `MappedNeuroVec`, or
  `FileBackedNeuroVec`).

- file:

  Optional output file name. If `NULL`, a temporary file with extension
  `.nii` is created.

- ...:

  Additional arguments passed to methods (e.g. `data_type`,
  `overwrite`).

- data_type:

  Character string specifying the output data type for the NIfTI file.
  Should be one of: "BINARY", "UBYTE", "SHORT", "INT", "FLOAT",
  "DOUBLE". Default is "FLOAT".

- overwrite:

  Logical; if `TRUE`, overwrite an existing file at the specified path.
  Default is `FALSE`.

## Value

A
[`MappedNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/MappedNeuroVec-class.md)
(or `x` itself if already memory-mapped).

A
[`MappedNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/MappedNeuroVec-class.md)
(or `x` itself if it is already memory-mapped).
