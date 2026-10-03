# Profile every evaluated vignette chunk; no cache or evaluation skips.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
repo <- normalizePath(args[1])
out <- normalizePath(args[2], mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (nzchar(Sys.getenv('PROFILE_LIB'))) .libPaths(c(Sys.getenv('PROFILE_LIB'), .libPaths()))
stopifnot(as.character(packageVersion('neuroim2')) == read.dcf(file.path(repo, 'DESCRIPTION'))[1, 'Version'])
file.copy(file.path(repo, 'vignettes'), out, recursive = TRUE)
src <- file.path(out, 'vignettes')
writeLines(capture.output(sessionInfo()), file.path(out, 'session-info.txt'))
chunks <- data.frame(vignette=character(), chunk=character(), elapsed=numeric())
totals <- data.frame(vignette=character(), elapsed=numeric())
starts <- new.env(parent=emptyenv())
knitr::knit_hooks$set(profile_runtime=function(before, options, envir) {
  if (before) {
    starts[[options$label]] <- proc.time()[['elapsed']]
  } else {
    seconds <- proc.time()[['elapsed']] - starts[[options$label]]
    chunks <<- rbind(chunks, data.frame(vignette=current, chunk=options$label, elapsed=seconds))
    write.csv(chunks, file.path(out, 'chunk-timings.csv'), row.names=FALSE)
    cat(sprintf('TIMING %s/%s %.3f seconds\n', current, options$label, seconds), file=stderr())
  }
  ''
})
knitr::opts_chunk$set(profile_runtime=TRUE, cache=FALSE)
for (input in list.files(src, '[.]Rmd$', full.names=TRUE)) {
  current <- basename(input)
  start <- proc.time()[['elapsed']]
  cat('RENDER', current, 'with', find.package('neuroim2'), '\n')
  Rprof(file.path(out, paste0(current, '.Rprof')), interval=0.02)
  rmarkdown::render(input, quiet=TRUE, envir=new.env(parent=globalenv()))
  Rprof(NULL)
  totals <- rbind(totals, data.frame(vignette=current, elapsed=proc.time()[['elapsed']]-start))
  write.csv(totals, file.path(out, 'vignette-timings.csv'), row.names=FALSE)
  prof <- summaryRprof(file.path(out, paste0(current, '.Rprof')))
  write.csv(head(prof$by.total,40), file.path(out, paste0(current,'.profile-total.csv')))
  write.csv(head(prof$by.self,40), file.path(out, paste0(current,'.profile-self.csv')))
}
print(totals[order(totals$elapsed, decreasing=TRUE), ])
print(head(chunks[order(chunks$elapsed, decreasing=TRUE), ], 20))
