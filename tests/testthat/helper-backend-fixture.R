# Independent native-endian NIfTI fixture; no package writer/reader in oracle.
backend_binary_fixture <- function(type = "SHORT", dims = c(3L, 2L, 1L, 3L)) {
  path <- tempfile(fileext = ".nii")
  stored <- switch(type,
    FLOAT = rep_len(c(-1.25, 0, 2.5, NA_real_, 4, 8), prod(dims)),
    SHORT = as.integer(seq_len(prod(dims)) - 9L),
    UBYTE = as.integer((seq_len(prod(dims)) * 7) %% 256))
  code <- switch(type, FLOAT = 16L, SHORT = 4L, UBYTE = 2L)
  bytes <- switch(type, FLOAT = 4L, SHORT = 2L, UBYTE = 1L)
  slope <- if (type == "SHORT") -2.5 else 1
  intercept <- if (type == "SHORT") 7 else 0
  con <- file(path, "w+b")
  on.exit(close(con))
  writeBin(raw(352), con)
  put <- function(offset, value, size) {
    seek(con, offset, origin = "start", rw = "write")
    writeBin(value, con, size = size, endian = .Platform$endian)
  }
  put(0, 348L, 4)
  put(40, c(4L, as.integer(dims), 1L, 1L, 1L), 2)
  put(70, code, 2)
  put(72, bytes * 8L, 2)
  put(76, rep(1, 8), 4)
  put(108, c(352, slope, intercept), 4)
  seek(con, 344, rw = "write")
  writeBin(as.raw(c(110, 43, 49, 0)), con)
  if (type == "UBYTE") stored <- as.raw(stored)
  put(352, stored, bytes)
  list(path = path, expected = array(as.numeric(stored) * slope + intercept, dims))
}
