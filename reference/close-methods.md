# Close a BinaryReader or BinaryWriter

Closes the underlying connection associated with a BinaryReader or
BinaryWriter object. This should be called when you're done with the
reader/writer to free system resources.

## Usage

``` r
# S4 method for class 'BinaryReader'
close(con)

# S4 method for class 'BinaryWriter'
close(con)
```

## Arguments

- con:

  The BinaryReader or BinaryWriter object to close.

## Value

Invisibly returns `NULL`, called for its side effect of closing the
connection.

## Examples

``` r
# \donttest{
# Create a temporary file and write some data
tmp <- tempfile()
writer <- BinaryWriter(tmp, byte_offset = 0L,
                      data_type = "DOUBLE", bytes_per_element = 8L)
write_elements(writer, rnorm(100))
close(writer)

# Read the data back
reader <- BinaryReader(tmp, byte_offset = 0L,
                      data_type = "DOUBLE", bytes_per_element = 8L)
data <- read_elements(reader, 100)
close(reader)

# Clean up
unlink(tmp)
# }
```
