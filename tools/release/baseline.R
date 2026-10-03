out <- Sys.getenv("RELEASE_OUT")
stopifnot(nzchar(out))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
baseline_sha <- "f1f7f00264b0566beb73d17654ff1e874e66bf3c"
src <- tempfile("master-source-")
lib <- tempfile("master-library-")
dir.create(src); dir.create(lib)
archive <- tempfile(fileext = ".tar")
stopifnot(system2("git", c("archive", baseline_sha, "-o", shQuote(archive))) == 0L)
untar(archive, exdir = src)
install.packages(src, lib = lib, repos = NULL, type = "source")
writeLines(baseline_sha, file.path(out, "BASELINE_SHA"))
for (mode in c("before", "after")) {
  args <- if (mode == "before") c(mode, lib) else mode
  log <- file.path(out, paste0("backend-regressions-", mode, ".log"))
  callr::rscript("tools/release/reproduce-backend-defects.R", cmdargs = args,
                 stdout = log, stderr = "2>&1")
  cat(readLines(log), sep = "\n")
}
