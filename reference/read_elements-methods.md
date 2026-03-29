# Read a sequence of elements from an input source

Read a specified number of elements from a
[BinaryReader](https://bbuchsbaum.github.io/neuroim2/reference/BinaryReader-class.md)
object.

## Usage

``` r
read_elements(x, num_elements)

# S4 method for class 'BinaryReader,numeric'
read_elements(x, num_elements)
```

## Arguments

- x:

  Object of class
  [BinaryReader](https://bbuchsbaum.github.io/neuroim2/reference/BinaryReader-class.md)

- num_elements:

  Integer specifying number of elements to read

## Value

A `vector` containing the elements read from `x`.

Numeric vector of read elements

## Examples

``` r
# Create a temporary binary file with test data
tmp <- tempfile()
con <- file(tmp, "wb")
test_data <- rnorm(100)
writeBin(test_data, con, size = 8)
close(con)

# Create a BinaryReader and read the data
reader <- BinaryReader(tmp, byte_offset = 0L,
                      data_type = "double", bytes_per_element = 8L)
data <- read_elements(reader, 100)
close(reader)

# Clean up
unlink(tmp)
# \donttest{
# Create a temporary binary file with some test data
tmp <- tempfile()
con <- file(tmp, "wb")
test_data <- rnorm(100)
writeBin(test_data, con, size = 8)
close(con)

# Create reader and read the data
reader <- BinaryReader(tmp, byte_offset = 0L,
                      data_type = "DOUBLE", bytes_per_element = 8L)
data <- read_elements(reader, 100)
close(reader)

# Clean up
unlink(tmp)
# }
```
