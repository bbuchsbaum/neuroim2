# Create Column Reader Object

Create a new instance of the
[ColumnReader](https://bbuchsbaum.github.io/neuroim2/reference/ColumnReader-class.md)
class for reading column-oriented data.

## Usage

``` r
ColumnReader(nrow, ncol, reader)
```

## Arguments

- nrow:

  Integer specifying number of rows in data

- ncol:

  Integer specifying number of columns in data

- reader:

  Function that takes column indices and returns matrix

## Value

An object of class
[ColumnReader](https://bbuchsbaum.github.io/neuroim2/reference/ColumnReader-class.md)

## Examples

``` r

reader_func <- function(cols) {
  matrix(rnorm(100 * length(cols)), 100, length(cols))
}
col_reader <- ColumnReader(nrow = 100L, ncol = 10L, reader = reader_func)
```
