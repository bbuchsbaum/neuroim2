# Run from the package root: Rscript --vanilla dev/issue-triage/probe.R
# Optional baseline: pass a Git revision, e.g. 77b1ddb.
pkgload::load_all(quiet = TRUE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  # Evaluate historical functions without changing the checkout or installed package.
  for (file in c("sparse_neurovec.R", "plot-overlay.R")) {
    src <- system2("git", c("show", paste0(args[1], ":R/", file)), stdout = TRUE)
    eval(parse(text = src), envir = .GlobalEnv)
  }
}
sp <- NeuroSpace(c(2, 2, 2, 3))
m <- array(seq_len(8) %in% c(1, 3, 6), c(2, 2, 2))
y <- matrix(seq_len(9), 3, 3)
v <- SparseNeuroVec(y, sp, m)
cat("#31 auto square is voxels x time:", isTRUE(all.equal(unname(series(v, c(1L,3L,6L))), unname(t(y)))), "\n")
cat("#5 impossible time rejected:", inherits(try(SparseNeuroVec(y, NeuroSpace(c(2,2,2,1000000)), m), silent=TRUE), "try-error"), "\n")
# #1: the retired MNI_SPACE_1MM dataset is no longer shipped. Exercise the
# documented world-to-grid workflow with a synthetic template, then render it.
s3 <- NeuroSpace(c(20, 20, 20), origin = c(-10, -10, -10))
bg <- NeuroVol(array(0, c(20,20,20)), s3)
roi_list <- list(c(0,0,0), c(3,3,3))
rois <- purrr::map(roi_list, ~spherical_roi(bg, coord_to_grid(s3, .x), radius=3, fill=1))
stopifnot(length(rois)==2, all(vapply(rois, function(x) length(x)>0, logical(1))))
cat("#1 world-coordinate ROI mapping succeeds; legacy dataset available:", "MNI_SPACE_1MM" %in% data(package="neuroim2")$results[,"Item"], "\n")
# Synthetic skull-stripped background with a small bright tissue component.
a <- array(0, c(20,20,1)); a[5:15,5:15,1] <- 100; a[9:10,9:10,1] <- 200
s <- NeuroSpace(dim(a)); bg <- NeuroVol(a,s); ov <- NeuroVol(array(0,dim(a)),s)
p <- plot_overlay(bg,ov,zlevels=1,draw=FALSE,assemble=FALSE)
cat("#19 default background limits:", p[[1]]$scales$get_scales("fill")$limits, "\n")
print(sessionInfo())
