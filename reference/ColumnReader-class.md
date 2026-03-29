# ColumnReader

A class that supports reading of data from a matrix-like storage format,
such as a file or a database, in a column-wise manner.

## Slots

- `nrow`:

  An `integer` representing the number of rows in the matrix-like
  storage.

- `ncol`:

  An `integer` representing the number of columns in the matrix-like
  storage.

- `reader`:

  A `function` that takes a set of column indices as input and returns a
  `matrix` containing the requested columns from the storage.
