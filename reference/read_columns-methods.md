# Read a set of column vector from an input source (e.g. `ColumnReader`)

Read a set of column vector from an input source (e.g. `ColumnReader`)

## Usage

``` r
read_columns(x, column_indices)

# S4 method for class 'ColumnReader,integer'
read_columns(x, column_indices)
```

## Arguments

- x:

  the input channel

- column_indices:

  the column indices

## Value

A numeric `matrix` consisting of the requested column vectors.

## Examples

``` r
# Create a reader function that returns random data
reader_func <- function(cols) {
  matrix(rnorm(100 * length(cols)), 100, length(cols))
}

# Create a ColumnReader with 100 rows and 10 columns
col_reader <- ColumnReader(nrow = 100L, ncol = 10L, reader = reader_func)

# Read columns 1, 3, and 5
cols <- read_columns(col_reader, c(1L, 3L, 5L))
dim(cols) == c(100, 3)
#> [1] TRUE TRUE
```
