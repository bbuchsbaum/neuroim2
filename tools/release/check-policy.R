# A successful R CMD check exit status is not a release decision: CRAN's
# pretest also enforces overall runtime. Keep this independent of rcmdcheck
# so the same policy can be applied to official external-service results.
review_check <- function(errors, warnings, notes, elapsed, service,
                         approvals = "tools/release/reviewed-notes.dcf") {
  stopifnot(length(service) == 1L, nzchar(service))
  limit <- 600
  normalize <- function(x) {
    # R CMD check can decorate a stage heading with CPU/elapsed timings.
    # Strip only that heading annotation; the NOTE body remains exact and
    # complete-check elapsed time is reviewed independently below.
    x <- sub("^(checking [^\\n]+ \\.\\.\\.) \\[[0-9.]+[smh]/[0-9.]+[smh]\\] NOTE(\\n|$)",
             "\\1 NOTE\\2", x, perl = TRUE)
    unname(gsub("[[:space:]]+", " ", trimws(x)))
  }
  approved <- read.dcf(approvals)
  note_ids <- character()
  unreviewed <- character()
  for (note in notes) {
    matching <- which(vapply(seq_len(nrow(approved)), function(i) {
      service %in% strsplit(approved[i, "Services"], " +")[[1]] &&
        identical(normalize(note), normalize(approved[i, "Note"]))
    }, logical(1)))
    if (length(matching) == 1L) {
      note_ids <- c(note_ids, approved[matching, "Id"])
    } else {
      unreviewed <- c(unreviewed, note)
    }
  }
  reasons <- c(
    if (length(errors)) sprintf("%d check error(s)", length(errors)),
    if (length(warnings)) sprintf("%d check warning(s)", length(warnings)),
    if (length(unreviewed)) sprintf("%d unreviewed NOTE(s)", length(unreviewed)),
    if (length(elapsed) != 1L || !is.finite(elapsed) || elapsed <= 0) {
      "Missing or invalid complete-check elapsed time"
    } else if (elapsed > limit) {
      sprintf("Overall checktime %.1f seconds exceeds %d seconds", elapsed, limit)
    }
  )
  list(pass = length(reasons) == 0L, reasons = reasons, note_ids = note_ids,
       unreviewed_notes = unreviewed, elapsed = elapsed, limit = limit,
       service = service)
}

write_check_policy <- function(review, result, out, source_sha, checksum,
                               log_file = "check/neuroim2.Rcheck/00check.log") {
  stopifnot(file.exists(file.path(out, log_file)))
  report <- c(
    Service = review$service, SourceSHA = source_sha, SHA256 = checksum,
    CheckedAtUTC = format(Sys.time(), tz = "UTC", usetz = TRUE),
    ElapsedSeconds = as.character(review$elapsed),
    LimitSeconds = as.character(review$limit),
    Errors = as.character(length(result$errors)),
    Warnings = as.character(length(result$warnings)),
    Notes = as.character(length(result$notes)),
    ReviewedNoteIDs = paste(review$note_ids, collapse = " "),
    Gate = if (review$pass) "PASS" else "FAIL",
    Reasons = paste(review$reasons, collapse = "; "), LogFile = log_file,
    LogSHA256 = digest::digest(file = file.path(out, log_file), algo = "sha256")
  )
  write.dcf(as.data.frame(as.list(report)), file.path(out, "check-policy.dcf"))
  saveRDS(review, file.path(out, "check-policy.rds"))
  if (!review$pass) stop(paste(review$reasons, collapse = "; "), call. = FALSE)
  invisible(report)
}
