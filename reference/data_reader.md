# Create a Data Reader

Creates a data reader for accessing neuroimaging data from various file
formats. The reader provides a unified interface for reading data
regardless of the underlying format.

## Usage

``` r
data_reader(x, offset)
```

## Arguments

- x:

  An object containing metadata required to create the reader (e.g.,
  file path, format info)

- offset:

  Numeric. Byte offset where data reading should begin. Default is 0.

## Value

A BinaryReader object configured for the specific data format

## Details

Create a Data Reader for Neuroimaging Data

The data_reader function is a generic that creates appropriate readers
for different neuroimaging formats. It handles:

- File format detection and validation

- Endianness configuration

- Data type conversion

- Compression handling (e.g., gzip)

- Proper byte alignment

## See also

[`read_header`](https://bbuchsbaum.github.io/neuroim2/reference/read_header.md)
for reading headers,
[`BinaryReader`](https://bbuchsbaum.github.io/neuroim2/reference/BinaryReader-class.md)
for reading binary data

## Examples

``` r

# Create reader for NIFTI file
meta <- read_header(system.file("extdata", "global_mask_v4.nii", package="neuroim2"))
reader <- data_reader(meta, offset = 0)

# Read first 100 voxels
data <- read_elements(reader, 100)
close(reader)

```
