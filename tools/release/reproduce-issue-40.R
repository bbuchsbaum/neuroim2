# Exact public-writer reproduction from issue #40, with the SHORT defect
# assertion replaced by dense/mapped agreement. Uses only selected indices.
library(neuroim2)

check_storage_type <- function(dtype) {
  expected <- array(seq_len(72) / 8 - 1.3, c(3L, 2L, 2L, 6L))
  path <- tempfile(fileext = ".nii")
  on.exit(unlink(path), add = TRUE)
  write_vec(DenseNeuroVec(expected, NeuroSpace(dim(expected))),
            path, data_type = dtype)
  header <- read_header(path)
  dense <- read_vec(path)
  mapped <- MappedNeuroVec(path)
  on.exit(mmap::munmap(mapped@filemap), add = TRUE, after = FALSE)
  idx <- as.double(c(1, 12, 13, 37, 72))
  dv <- linear_access(dense, idx)
  mv <- linear_access(mapped, idx)
  cat(dtype, "slope =", header@slope, "intercept =", header@intercept, "\n")
  print(data.frame(index = idx, expected = as.vector(expected)[idx],
                   dense = dv, mapped = mv))
  cat("max dense/mapped difference:", max(abs(dv - mv)), "\n")
  stopifnot(max(abs(dv - as.vector(expected)[idx])) < 1e-3)
  stopifnot(max(abs(mv - as.vector(expected)[idx])) < 1e-3)
  if (dtype == "FLOAT") stopifnot(identical(dv, mv))
  if (dtype == "SHORT") stopifnot(isTRUE(all.equal(dv, mv, tolerance = 1e-12)))
}
check_storage_type("FLOAT")
check_storage_type("SHORT")
print(packageVersion("neuroim2"))
print(sessionInfo())
