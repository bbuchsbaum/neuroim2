# Run from the repository root. The output must be outside the source tree.
out <- Sys.getenv("RELEASE_OUT", unset = "")
stopifnot(nzchar(out))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
Sys.setenv(NOT_CRAN = "false")
writeLines(system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
           file.path(out, "SOURCE_SHA"))
writeLines(capture.output(sessionInfo()), file.path(out, "session-info.txt"))

# Both options are explicit: no --no-manual or --no-build-vignettes shortcuts.
tarball <- pkgbuild::build(".", dest_path = out, vignettes = TRUE,
                          manual = TRUE, quiet = FALSE)
writeLines(system2("sha256sum", shQuote(tarball), stdout = TRUE),
           file.path(out, "SHA256SUMS"))
result <- rcmdcheck::rcmdcheck(tarball, args = "--as-cran",
                             check_dir = file.path(out, "check"),
                             error_on = "never")
saveRDS(result[c("errors", "warnings", "notes")], file.path(out, "results.rds"))
writeLines(capture.output(print(result)), file.path(out, "results.txt"))
stopifnot(length(result$errors) == 0L, length(result$warnings) == 0L)
