# Build with current release/patched R, per CRAN submission policy.
# Run from the repository root. The output must be outside the source tree.
stopifnot(!grepl("development", R.version$status, ignore.case = TRUE))
out <- Sys.getenv("RELEASE_OUT", unset = "")
stopifnot(nzchar(out))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
Sys.setenv(NOT_CRAN = "false")
writeLines(system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
           file.path(out, "SOURCE_SHA"))
writeLines(capture.output(sessionInfo()), file.path(out, "build-session-info.txt"))

# Both options are explicit: no --no-manual or --no-build-vignettes shortcuts.
tarball <- pkgbuild::build(".", dest_path = out, vignettes = TRUE,
                          manual = TRUE, quiet = FALSE)
checksum <- if (nzchar(Sys.which("sha256sum"))) {
  system2("sha256sum", shQuote(tarball), stdout = TRUE)
} else {
  system2("shasum", c("-a", "256", shQuote(tarball)), stdout = TRUE)
}
stopifnot(is.null(attr(checksum, "status")), length(checksum) == 1L)
writeLines(checksum, file.path(out, "SHA256SUMS"))
source("tools/release/check-vignette-assets.R")
unpacked <- file.path(out, "source-docs")
dir.create(unpacked)
untar(tarball, files = grep("^neuroim2/inst/doc/", untar(tarball, list = TRUE),
                           value = TRUE), exdir = unpacked)
source_assets <- check_vignette_assets(file.path(unpacked, "neuroim2/inst/doc"))
write.csv(source_assets, file.path(out, "built-vignette-assets.csv"), row.names = FALSE)
