# BinaryReader Class

Class supporting reading of bulk binary data from a connection

## Slots

- `input`:

  The binary input connection

- `byte_offset`:

  The number of bytes to skip at the start of input

- `data_type`:

  The data type of the binary elements

- `bytes_per_element`:

  The number of bytes in each data element (e.g. 4 or 8 for floating
  point numbers)

- `endian`:

  The endianness of the binary input connection

- `signed`:

  Logical indicating whether the data type is signed
