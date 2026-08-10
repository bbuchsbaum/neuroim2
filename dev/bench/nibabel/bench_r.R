suppressMessages(library(neuroim2))

tm <- function(label, expr, reps = 3) {
  f <- function() { force(expr) }
  # warmup
  invisible(eval.parent(substitute(expr)))
  ts <- numeric(reps)
  for (i in seq_len(reps)) {
    t0 <- proc.time()[["elapsed"]]
    invisible(eval.parent(substitute(expr)))
    ts[i] <- proc.time()[["elapsed"]] - t0
  }
  cat(sprintf("%-46s %8.3f s\n", label, median(ts)))
  invisible(median(ts))
}

B <- Sys.getenv("NEUROIM2_BENCH_DIR", "/tmp/bench")

tm("header  3D .nii (read_header)",        read_header(file.path(B, "vol3d_i16.nii")))
tm("header  4D .nii (read_header)",        read_header(file.path(B, "fmri_i16.nii")))
tm("read    3D int16 .nii (read_vol)",     read_vol(file.path(B, "vol3d_i16.nii")))
tm("read    3D float32 .nii (read_vol)",   read_vol(file.path(B, "vol3d_f32.nii")))
tm("read    3D float32 .nii.gz (read_vol)",read_vol(file.path(B, "vol3d_f32.nii.gz")))
tm("read    4D int16 .nii (read_vec)",     read_vec(file.path(B, "fmri_i16.nii")), reps = 3)
tm("read    4D int16 .nii.gz (read_vec)",  read_vec(file.path(B, "fmri_i16.nii.gz")), reps = 2)
tm("read    4D subvolume k=100 (read_vol)",read_vol(file.path(B, "fmri_i16.nii"), index = 100))

v <- read_vol(file.path(B, "vol3d_f32.nii"))
tm("write   3D float32 .nii (write_vol)",  write_vol(v, tempfile(fileext = ".nii")))
tm("write   3D float32 .nii.gz",           write_vol(v, tempfile(fileext = ".nii.gz")), reps = 2)

x <- read_vec(file.path(B, "fmri_i16.nii"))
idx <- as.integer(seq(1, prod(dim(x)[1:3]), length.out = 1000))
tm("series  1000 voxels x 200 tp (in-mem)", series(x, idx), reps = 5)

cat("\n-- lazy/mmap 4D access --\n")
ok <- tryCatch({
  xm <- read_vec(file.path(B, "fmri_i16.nii"), mode = "mmap")
  tm("mmap    open 4D", read_vec(file.path(B, "fmri_i16.nii"), mode = "mmap"))
  tm("mmap    series 1000 voxels x 200 tp", series(xm, idx), reps = 3)
  TRUE
}, error = function(e) { cat("mmap mode failed:", conditionMessage(e), "\n"); FALSE })
