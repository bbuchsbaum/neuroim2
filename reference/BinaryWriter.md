# Create Binary Writer Object

Create a new instance of the
[BinaryWriter](https://bbuchsbaum.github.io/neuroim2/reference/BinaryWriter-class.md)
class for writing bulk binary data.

## Usage

``` r
BinaryWriter(
  output,
  byte_offset,
  data_type,
  bytes_per_element,
  endian = .Platform$endian
)
```

## Arguments

- output:

  Character string (file name) or connection object to write to

- byte_offset:

  Integer specifying bytes to skip at start of output

- data_type:

  Character string specifying R data type ('integer', 'double', etc.)

- bytes_per_element:

  Integer specifying bytes per data element (e.g., 4 or 8)

- endian:

  Character string specifying endianness ('big' or 'little', default:
  platform-specific)

## Value

An object of class
[BinaryWriter](https://bbuchsbaum.github.io/neuroim2/reference/BinaryWriter-class.md)

## See also

[`BinaryReader`](https://bbuchsbaum.github.io/neuroim2/reference/BinaryReader.md)
for reading binary data

## Examples

``` r
# \donttest{

tmp <- tempfile()
# Write to existing connection with offset
con <- file(tmp, "wb")
writer <- BinaryWriter(con, byte_offset = 100L,
                      data_type = "integer", bytes_per_element = 4L)
close(writer)
unlink(tmp)
# }
```
