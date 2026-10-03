# Rscript tools/release/final-gate.R ARTIFACT_DIR CHECK_DIR [CHECK_DIR ...]
# Check directories must contain machine-generated check-policy.dcf and the
# unchanged log it identifies. This command never uploads or submits anything.
args <- commandArgs(TRUE)
stopifnot(length(args) >= 2L)
candidate <- normalizePath(args[1], mustWork=TRUE)
tarball <- list.files(candidate, "^neuroim2_.*[.]tar[.]gz$", full.names=TRUE)
stopifnot(length(tarball) == 1L)
sha <- digest::digest(file=tarball, algo="sha256")
source_sha <- readLines(file.path(candidate, "SOURCE_SHA"))
recorded <- strsplit(readLines(file.path(candidate, "SHA256SUMS")), " +")[[1]][1]
stopifnot(identical(sha, recorded), length(source_sha) == 1L,
          grepl("^[0-9a-f]{40}$", source_sha))
required <- c("github-r-devel", "winbuilder-r-devel", "rhub-ubuntu-clang",
              "rhub-clang-ubsan")
reports <- lapply(args[-1], function(dir) {
  x <- read.dcf(file.path(dir, "check-policy.dcf"))
  stopifnot(nrow(x) == 1L, x[1,"SHA256"] == sha,
            x[1,"SourceSHA"] == source_sha, x[1,"Gate"] == "PASS",
            x[1,"Errors"] == "0", x[1,"Warnings"] == "0",
            as.numeric(x[1,"ElapsedSeconds"]) > 0,
            as.numeric(x[1,"ElapsedSeconds"]) <= 600,
            as.numeric(x[1,"LimitSeconds"]) == 600)
  log <- normalizePath(file.path(dir, x[1,"LogFile"]), mustWork=TRUE)
  stopifnot(digest::digest(file=log, algo="sha256") == x[1,"LogSHA256"])
  x
})
services <- vapply(reports, function(x) x[1,"Service"], "")
stopifnot(!anyDuplicated(services), all(required %in% services))
evidence <- do.call(rbind, reports)
write.dcf(evidence, file.path(candidate, "final-gate-evidence.dcf"))
cat("PASS: complete checks within ten minutes, reviewed notes, matching source",
    "and tarball on all required services.\nSHA256:", sha, "\n")
