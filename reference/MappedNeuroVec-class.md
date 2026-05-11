# MappedNeuroVec Class

A class representing a four-dimensional brain image backed by a
memory-mapped file. This class provides efficient access to large brain
images without loading the entire dataset into memory.

The `MappedNeuroVec` class provides memory-efficient access to large
neuroimaging datasets through memory mapping. This allows processing of
datasets larger than available RAM by keeping data on disk and only
loading requested portions into memory.

Creates a `MappedNeuroVec` object that provides efficient, memory-mapped
access to large neuroimaging datasets. This allows processing of data
larger than available RAM by keeping it on disk and only loading
requested portions into memory.

## Usage

``` r
MappedNeuroVec(file_name, label = basename(file_name))
```

## Arguments

- file_name:

  Character string specifying the path to the neuroimaging file.
  Supported formats include NIFTI (.nii) and ANALYZE (.hdr/.img).

- label:

  Optional character string providing a label for the vector

## Value

A new `MappedNeuroVec` object providing:

- Memory-mapped access to the data

- Spatial and temporal indexing

- Efficient data extraction

- Automatic memory management

## Details

MappedNeuroVec objects use memory-mapped files to store and access large
4D brain images efficiently. This approach allows for rapid access to
specific portions of the data without requiring the entire dataset to be
loaded into memory at once.

Create a Memory-Mapped Neuroimaging Vector

The function implements several key features:

- Zero-copy access to file data

- Automatic memory management

- Support for large datasets

- Efficient random access

- Proper cleanup on object deletion

Memory mapping is particularly useful when:

- Working with large datasets

- Only portions of data are needed at once

- Random access is required

- Multiple processes need to share data

## Slots

- `filemap`:

  An object of class `mmap` representing the memory-mapped file
  containing the brain image data.

- `offset`:

  An integer representing the byte offset within the memory-mapped file
  where the brain image data starts.

## Methods

This class inherits methods from
[`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md)
and implements the `ArrayLike4D` interface. Additional methods specific
to memory-mapped operations may be available.

## Implementation Details

The class uses the `mmap` package to establish a memory mapping between
the file and memory space. Key features include:

- Zero-copy access to file data

- Automatic memory management

- Support for large datasets

- Efficient random access

## See also

[`NeuroVec-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.md)
for the parent class. [`mmap`](https://rdrr.io/pkg/mmap/man/mmap.html)
for details on memory-mapped file objects.

`MappedNeuroVec` for creating instances of this class

[`mmap`](https://rdrr.io/pkg/mmap/man/mmap.html) for memory mapping
details

## Examples

``` r

# Create a MappedNeuroVec object (pseudo-code)
file_path <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
mapped_vec <- MappedNeuroVec(file_path)

# Access a subset of the data
subset <- mapped_vec[,,, 1:2]


# Create mapped vector from NIFTI file
mvec <- MappedNeuroVec(system.file("extdata", "global_mask_v4.nii", package = "neuroim2"))

# Extract first volume
vol1 <- mvec[[1]]

# Get dimensions
dim(mvec)
#> [1] 64 64 25  4

# Access specific timepoint
timepoint <- mvec[, , , 2]

```
