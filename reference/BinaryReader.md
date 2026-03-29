# Create Binary Reader Object

Create a new instance of the
[BinaryReader](https://bbuchsbaum.github.io/neuroim2/reference/BinaryReader-class.md)
class for reading bulk binary data.

## Usage

``` r
BinaryReader(
  input,
  byte_offset,
  data_type,
  bytes_per_element,
  endian = .Platform$endian,
  signed = TRUE
)
```

## Arguments

- input:

  Character string (file name) or connection object to read from

- byte_offset:

  Integer specifying bytes to skip at start of input

- data_type:

  Character string specifying R data type ('integer', 'double', etc.)

- bytes_per_element:

  Integer specifying bytes per data element (e.g., 4 or 8)

- endian:

  Character string specifying endianness ('big' or 'little', default:
  platform-specific)

- signed:

  Logical indicating if data type is signed (default: TRUE)

## Value

An object of class
[BinaryReader](https://bbuchsbaum.github.io/neuroim2/reference/BinaryReader-class.md)

## See also

[`BinaryWriter`](https://bbuchsbaum.github.io/neuroim2/reference/BinaryWriter.md)
for writing binary data

## Examples

``` r
# \donttest{
# Create a temporary binary file
tmp <- tempfile()
writeBin(rnorm(100), tmp, size = 8)


# Read from existing connection with offset
con <- file(tmp, "rb")
reader <- BinaryReader(con, byte_offset=0,
                      data_type = "DOUBLE", bytes_per_element = 8L)
close(reader)

# Clean up
unlink(tmp)
# }
```
