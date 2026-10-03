# Check the exact release-R-built source tarball with current R-devel.
out <- Sys.getenv("RELEASE_OUT", unset = "")
stopifnot(nzchar(out))
out <- normalizePath(out, mustWork = TRUE)
Sys.setenv(NOT_CRAN = "false")
head <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
stopifnot(identical(readLines(file.path(out, "SOURCE_SHA")), head))
writeLines(capture.output(sessionInfo()), file.path(out, "check-session-info.txt"))
tarball <- list.files(out, "^neuroim2_.*[.]tar[.]gz$", full.names = TRUE)
stopifnot(length(tarball) == 1L)
checksum <- if (nzchar(Sys.which("sha256sum"))) {
  system2("sha256sum", shQuote(tarball), stdout = TRUE)
} else {
  system2("shasum", c("-a", "256", shQuote(tarball)), stdout = TRUE)
}
stopifnot(is.null(attr(checksum, "status")), length(checksum) == 1L)
recorded <- readLines(file.path(out, "SHA256SUMS"))
stopifnot(identical(strsplit(checksum, " ")[[1]][1],
                    strsplit(recorded, " ")[[1]][1]))
source("tools/release/check-vignette-assets.R")
source_assets <- check_vignette_assets(file.path(out, "source-docs/neuroim2/inst/doc"))
started <- proc.time()[["elapsed"]]
result <- rcmdcheck::rcmdcheck(tarball, args = c("--as-cran", "--timings"),
                             check_dir = file.path(out, "check"),
                             error_on = "never")
elapsed <- proc.time()[["elapsed"]] - started
saveRDS(result[c("errors", "warnings", "notes")], file.path(out, "results.rds"))
writeLines(capture.output(print(result)), file.path(out, "results.txt"))
source("tools/release/check-policy.R")
review <- review_check(result$errors, result$warnings, result$notes, elapsed,
                      Sys.getenv("RELEASE_CHECK_SERVICE", "local"))
if (!review$pass) write_check_policy(review, result, out, head,
                                     strsplit(recorded, " ")[[1]][1])
installed_assets <- check_vignette_assets(
  file.path(out, "check/neuroim2.Rcheck/neuroim2/doc"))
write.csv(rbind(source_assets, installed_assets),
          file.path(out, "vignette-assets.csv"), row.names = FALSE)
print(rbind(source_assets, installed_assets))
write_check_policy(review, result, out, head,
                   strsplit(recorded, " ")[[1]][1])
