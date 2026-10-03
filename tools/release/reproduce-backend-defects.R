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

mask <- array(FALSE, c(4, 4, 4)); mask[1:4] <- TRUE
x <- SparseNeuroVec(rbind(c(1, 3, 5, 7), c(10, 30, 50, 70)),
                    NeuroSpace(c(4, 4, 4, 2)), mask,
                    orientation = "time_x_voxels")
y <- downsample(x, factor = 0.5)
actual <- series(y, indices(y))
expected <- rbind(c(2, 6), c(20, 60))
cat("Downsample expected:\n"); print(expected)
cat("Downsample observed:\n"); print(actual)
downsample_ok <- isTRUE(all.equal(actual, expected))

mask <- array(FALSE, c(2, 2, 2)); mask[c(1, 3, 6)] <- TRUE
vox <- c(1L, 3L, 6L)
input <- matrix(c(1, 3, 9, 10, 20, 50, 11, 30, 80), 3, 3)
make <- function(m) SparseNeuroVec(m, NeuroSpace(c(2, 2, 2, nrow(m))),
                                  mask, orientation = "time_x_voxels")
x <- make(input)
scaling_ok <- isTRUE(all.equal(series(scale_series(x, FALSE, FALSE), vox), input))
arithmetic_ok <- isTRUE(all.equal(series(x + x, vox), input * 2))
concat_ok <- isTRUE(all.equal(series(concat(make(input[1, , drop = FALSE]),
                                          make(input[2:3, , drop = FALSE])), vox), input))
roi <- ROIVec(space(x), arrayInd(vox, c(2, 2, 2)), input)
roi_ok <- isTRUE(all.equal(series(as(roi, "SparseNeuroVec"), vox), input))
cat("Sequence correct:", sequence_ok, "\nMapped correct:", mapped_ok,
    "\nDownsample correct:", downsample_ok, "\nSparse scaling correct:", scaling_ok,
    "\nSparse arithmetic correct:", arithmetic_ok, "\nSparse concat correct:", concat_ok,
    "\nROI conversion correct:", roi_ok, "\n")
ok <- c(sequence_ok, mapped_ok, downsample_ok, scaling_ok, arithmetic_ok, concat_ok, roi_ok)
if (args[1] == "before") stopifnot(!any(ok))
if (args[1] == "after") stopifnot(all(ok))
