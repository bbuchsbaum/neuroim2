suppressMessages(library(neuroim2))
med <- function(f, reps=5) { f(); median(replicate(reps, system.time(f())[["elapsed"]])) }
set.seed(1)
D <- c(91L,109L,91L); SP <- c(2,2,2); R <- 8
g <- expand.grid(x=1:D[1],y=1:D[2],z=1:D[3])
mv <- ((g$x-46)/38)^2+((g$y-55)/48)^2+((g$z-46)/38)^2 <= 1
mask <- LogicalNeuroVol(array(mv, D), NeuroSpace(D, SP))
midx <- which(mv); K <- length(midx)
grid <- index_to_grid(mask, midx)
NC <- 1500L; set.seed(2); ctr <- grid[sort(sample(nrow(grid), NC)),,drop=FALSE]

t_roi <- med(function() for (i in seq_len(NC)) spherical_roi(mask, ctr[i,], R))
cat(sprintf("spherical_roi   %7.1f us/ROI   -> %6.1f s over %d voxels\n", 1e6*t_roi/NC, t_roi/NC*K, K))

# anisotropic
SPa <- c(1,1,4)
maska <- LogicalNeuroVol(array(mv, D), NeuroSpace(D, SPa))
ta <- med(function() for (i in seq_len(400L)) spherical_roi(maska, ctr[i,], R))
cat(sprintf("spherical_roi (1/1/4mm) %7.1f us/ROI\n", 1e6*ta/400))

# dense series
NT <- 200L; D4 <- c(48L,56L,40L,NT)
dv <- NeuroVec(array(rnorm(prod(D4)), D4), NeuroSpace(D4, c(3,3,3)))
roi <- coords(spherical_roi(NeuroSpace(D4[1:3], c(3,3,3)), c(24L,28L,20L), 8))
storage.mode(roi) <- "integer"
t_sd <- med(function() for (k in 1:200) series(dv, roi))/200*1e3
cat(sprintf("series(DenseNeuroVec, matrix)  %7.3f ms/call  (%d voxels x %d tp)\n", t_sd, nrow(roi), NT))

# sparse series
NT2 <- 250L
set.seed(5)
m2 <- array(mv, D)
smat <- matrix(rnorm(K*NT2), NT2, K)
sv <- SparseNeuroVec(smat, NeuroSpace(c(D,NT2), SP), LogicalNeuroVol(m2, NeuroSpace(D,SP)))
lins <- lapply(seq_len(500L), function(i) as.integer(grid_to_index(mask, coords(spherical_roi(mask, ctr[i,], R)))))
t_ss <- med(function() for (l in lins) series(sv, l))/500*1e6
cat(sprintf("series(SparseNeuroVec, idx)    %7.1f us/ROI  -> %5.1f s over %d voxels\n", t_ss, t_ss*K/1e6, K))
