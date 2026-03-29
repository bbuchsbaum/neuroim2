# Write a sequence of elements from an input source

Write a sequence of elements from an input source

## Usage

``` r
write_elements(x, els)

# S4 method for class 'BinaryWriter,numeric'
write_elements(x, els)
```

## Arguments

- x:

  the output channel

- els:

  the elements to write

## Value

Invisibly returns `NULL` after writing the elements.

## Examples

``` r
# Create a temporary binary file for writing
tmp <- tempfile()
writer <- BinaryWriter(tmp, byte_offset = 0L,
                      data_type = "DOUBLE", bytes_per_element = 8L)

# Write some random data
data <- rnorm(100)
write_elements(writer, data)
close(writer)

# Read back the data to verify
reader <- BinaryReader(tmp, byte_offset = 0L,
                      data_type = "double", bytes_per_element = 8L)
read_data <- read_elements(reader, 100)
close(reader)

# Verify data was written correctly
all.equal(data, read_data)
#> [1] TRUE

# Clean up
unlink(tmp)
# \donttest{
# Create a temporary binary file for writing
tmp <- tempfile()
writer <- BinaryWriter(tmp, byte_offset = 0L,
                      data_type = "DOUBLE", bytes_per_element = 8L)
# Write some data
write_elements(writer, rnorm(100))
close(writer)

# Clean up
unlink(tmp)
# }
```
