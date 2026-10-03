source("tools/release/check-policy.R")
check <- function(elapsed = 300, notes = character(), service = "github-r-devel",
                  errors = character(), warnings = character()) {
  review_check(errors, warnings, notes, elapsed, service)
}
stopifnot(check()$pass, check(600)$pass)
# Regression: archived CRAN pretest ended Status: OK but took 37 minutes.
stopifnot(!check(37 * 60)$pass, !check(600.01)$pass)
stopifnot(!check(NA_real_)$pass, !check(numeric())$pass,
          !check(-1)$pass, !check(Inf)$pass)
approved <- read.dcf("tools/release/reviewed-notes.dcf")[1, "Note"]
stopifnot(check(notes = approved)$pass)
# R-hub enables timings for every check stage. The same reviewed NOTE must
# match whether its heading is decorated, without discarding any NOTE body.
v8 <- read.dcf("tools/release/reviewed-notes.dcf")[3, "Note"]
timed_v8 <- sub("... NOTE", "... [2s/2s] NOTE", v8, fixed = TRUE)
stopifnot(check(notes = timed_v8, service = "rhub-ubuntu-clang")$pass,
          !check(notes = paste(timed_v8, "New diagnostic"))$pass,
          !check(notes = sub("unavailable", "failed validation", timed_v8))$pass,
          !check(notes = sub("HTML version of manual", "tests", timed_v8))$pass,
          !check(601, notes = timed_v8)$pass,
          !check(notes = timed_v8, service = "winbuilder-r-devel")$pass)
fortify <- read.dcf("tools/release/reviewed-notes.dcf")[4, "Note"]
stopifnot(check(notes = fortify, service = "rhub-clang-ubsan")$pass,
          !check(notes = fortify, service = "github-r-devel")$pass,
          !check(notes = fortify, service = "winbuilder-r-devel")$pass,
          !check(notes = sub("FORTIFY_SOURCE=3", "FORTIFY_SOURCE=2", fortify),
                 service = "rhub-clang-ubsan")$pass)
# An exception cannot hide a different/new note, timing problem, warning or
# a missing HTML validation step on the official Windows service.
stopifnot(!check(notes = paste(approved, "New diagnostic"))$pass,
          !check(notes = "checking tests ... NOTE\nExamples with CPU time > 5s")$pass,
          !check(notes = approved, service = "winbuilder-r-devel")$pass,
          !check(warnings = "warning")$pass, !check(errors = "error")$pass,
          !check(601, notes = approved)$pass)
cat("Release policy regression checks passed.\n")
