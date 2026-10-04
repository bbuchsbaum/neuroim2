# Rscript tools/release/record-winbuilder.R ARTIFACT_DIR RESULT_DIR
# RESULT_DIR contains downloaded 00check.log, upload.dcf and timing.dcf.
# upload.dcf records Service, SourceSHA, SHA256 and ReceiptURL at upload time.
# timing.dcf records StartUTC, EndUTC and Basis (e.g. check-log start and its
# HTTP Last-Modified, a conservative upper bound if server publication lags).
# This never uploads, resubmits, or communicates with CRAN.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
candidate <- normalizePath(args[1], mustWork=TRUE)
out <- normalizePath(args[2], mustWork=TRUE)
tarball <- list.files(candidate, "^neuroim2_.*[.]tar[.]gz$", full.names=TRUE)
stopifnot(length(tarball) == 1L)
sha <- digest::digest(file=tarball, algo="sha256")
source_sha <- readLines(file.path(candidate, "SOURCE_SHA"))
upload <- read.dcf(file.path(out, "upload.dcf"))
stopifnot(upload[1,"Service"] == "winbuilder-r-devel",
          upload[1,"SourceSHA"] == source_sha, upload[1,"SHA256"] == sha,
          grepl("^https://win-builder[.]r-project[.]org/", upload[1,"ReceiptURL"]))
log <- readLines(file.path(out,"00check.log"), warn=FALSE)
stopifnot(any(grepl("^\\* using R Under development",log)),
          any(grepl("^\\* using platform: x86_64-w64-mingw32",log)),
          any(log == "* DONE"), any(grepl("^Status:",log)))
required <- c("checking examples", "checking tests",
              "checking re-building of vignette outputs",
              "checking PDF version of manual", "checking HTML version of manual")
stopifnot(all(vapply(required, function(x) {
  lines <- grep(x, log, fixed=TRUE, value=TRUE)
  length(lines)==1L && !grepl("SKIPPED|NONE",lines)
},logical(1))))
timing <- read.dcf(file.path(out,"timing.dcf"))
start <- as.POSIXct(timing[1,"StartUTC"],format="%Y-%m-%d %H:%M:%S",tz="UTC")
end <- as.POSIXct(timing[1,"EndUTC"],format="%Y-%m-%d %H:%M:%S",tz="UTC")
stopifnot(nzchar(timing[1,"Basis"]),
          any(grepl(paste0("current time: ",format(start,tz="UTC")," UTC"),log,fixed=TRUE)))
elapsed <- as.numeric(difftime(end,start,units="secs"))
result <- rcmdcheck::parse_check(file=file.path(out,"00check.log"))
saveRDS(result[c("errors","warnings","notes")],file.path(out,"results.rds"))
source("tools/release/check-policy.R")
review <- review_check(result$errors,result$warnings,result$notes,elapsed,"winbuilder-r-devel")
write_check_policy(review,result,out,source_sha,sha,log_file="00check.log")
