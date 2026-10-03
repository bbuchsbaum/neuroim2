args <- commandArgs(TRUE)
stopifnot(length(args) >= 1L, args[1] %in% c("before", "after"))
if (length(args) > 1L) .libPaths(c(args[2], .libPaths()))
library(neuroim2)
cat("neuroim2:", as.character(packageVersion("neuroim2")), "\n")

# A non-symmetric square time-by-voxel result must keep its orientation.
x <- DenseNeuroVec(array(1:24, c(2, 2, 2, 3)), NeuroSpace(c(2, 2, 2, 3)))
y <- DenseNeuroVec(array(101:124, c(2, 2, 2, 3)), NeuroSpace(c(2, 2, 2, 3)))
expected <- rbind(series(x, 1:3), series(y, 1:3))
actual <- series(NeuroVecSeq(x, y), 1:3)
cat("Sequence expected:\n"); print(expected)
cat("Sequence observed:\n"); print(actual)
sequence_ok <- isTRUE(all.equal(actual, expected))

# Header and stored payload are independent of neuroim2's NIfTI writer.
source("tests/testthat/helper-backend-fixture.R")
fixture <- backend_binary_fixture("SHORT")
mapped <- MappedNeuroVec(fixture$path)
actual <- linear_access(mapped, 1:4)
expected <- as.vector(fixture$expected)[1:4]
cat("Mapped expected:", expected, "\nMapped observed:", actual, "\n")
mapped_ok <- isTRUE(all.equal(actual, expected))
mmap::munmap(mapped@filemap)
unlink(fixture$path)
cat("Sequence correct:", sequence_ok, "\nMapped correct:", mapped_ok, "\n")
if (args[1] == "before") stopifnot(!sequence_ok, !mapped_ok)
if (args[1] == "after") stopifnot(sequence_ok, mapped_ok)
