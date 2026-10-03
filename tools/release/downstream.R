# Isolated installs and subprocesses ensure each comparison uses its stated
# neuroim2 version, even after pak installs downstream dependencies.
out <- Sys.getenv("DOWNSTREAM_OUT")
stopifnot(nzchar(out))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
repos <- c(CRAN = "https://cloud.r-project.org")
options(repos = repos, timeout = 600)
Sys.setenv(NOT_CRAN = "false")
db <- available.packages(repos = repos, type = "source")
write.csv(db, file.path(out, "cran-package-index.csv"))
revdeps <- tools::package_dependencies("neuroim2", db, reverse = TRUE,
  which = c("Depends", "Imports", "LinkingTo", "Suggests"))[[1]]
writeLines(revdeps, file.path(out, "cran-reverse-dependencies.txt"))
stopifnot("bidser" %in% revdeps)
download <- function(pkg) download.packages(pkg, destdir = out, repos = repos,
                                           type = "source")[1, 2]
baseline <- download("neuroim2")
candidate <- list.files(Sys.getenv("RELEASE_OUT"), "^neuroim2_.*[.]tar[.]gz$",
                        full.names = TRUE)
stopifnot(length(candidate) == 1L)
libs <- file.path(out, c("baseline-lib", "candidate-lib"))
invisible(lapply(libs, dir.create))
for (i in 1:2) install.packages(c(baseline, candidate)[i], lib = libs[i],
                               repos = NULL, type = "source")
versions <- vapply(libs, function(lib) as.character(packageVersion("neuroim2", lib)), "")
writeLines(paste(c("baseline", "candidate"), versions), file.path(out, "versions.txt"))
stopifnot(versions[2] == read.dcf("DESCRIPTION")[1, "Version"])

new_failures <- character()
for (pkg in revdeps) {
  tarball <- download(pkg)
  src <- file.path(out, paste0(pkg, "-source"))
  dir.create(src)
  untar(tarball, exdir = src)
  # read.dcf's explicit fields retain absent optional entries as NA, as
  # package_dependencies requires every requested column to exist.
  desc <- read.dcf(file.path(src, pkg, "DESCRIPTION"),
    fields = c("Package", "Depends", "Imports", "LinkingTo", "Suggests"))
  deps <- tools::package_dependencies(pkg, db = desc,
    which = c("Depends", "Imports", "LinkingTo", "Suggests"))[[1]]
  pak::pkg_install(setdiff(deps, c("R", "neuroim2")), ask = FALSE)
  checks <- lapply(seq_along(libs), function(i) {
    rcmdcheck::rcmdcheck(tarball, args = c("--as-cran", "--no-manual"),
      check_dir = file.path(out, paste0(pkg, "-", c("baseline", "candidate")[i])),
      libpath = c(libs[i], .libPaths()), error_on = "never")
  })
  for (i in 1:2) writeLines(capture.output(print(checks[[i]])),
    file.path(out, paste0(pkg, "-", c("baseline", "candidate")[i], ".txt")))
  # Full logs retain existing findings. New errors/warnings fail this gate.
  before <- c(checks[[1]]$errors, checks[[1]]$warnings)
  after <- c(checks[[2]]$errors, checks[[2]]$warnings)
  new_failures <- c(new_failures, setdiff(after, before))
}

# A pinned, narrow integration check complements the complete CRAN checks.
# This does not claim a full fmridataset package check or all lab dependencies.
sha <- "ef141af6623b699168fa4329f304992c8b782ed0"
lab_tar <- file.path(out, "fmridataset-source.tar.gz")
download.file(paste0("https://api.github.com/repos/bbuchsbaum/fmridataset/tarball/", sha),
              lab_tar, mode = "wb")
lab_root <- file.path(out, "lab-source")
dir.create(lab_root)
untar(lab_tar, exdir = lab_root)
lab_src <- list.dirs(lab_root, recursive = FALSE)[1]
pak::local_install_deps(lab_src, dependencies = c("Depends", "Imports", "LinkingTo"), ask = FALSE)
writeLines(sha, file.path(out, "fmridataset-SOURCE_SHA"))
lab_result <- callr::r(function(src, lib, version, log) {
  .libPaths(c(lib, .libPaths()))
  stopifnot(as.character(packageVersion("neuroim2")) == version)
  sink(log, split = TRUE)
  on.exit(sink())
  testthat::test_local(src, filter = "nifti-array-source|feature-space",
                       reporter = "summary", stop_on_failure = TRUE)
  TRUE
}, args = list(lab_src, libs[2], versions[2], file.path(out, "fmridataset-tests.txt")),
   show = TRUE)
writeLines(new_failures, file.path(out, "new-cran-errors-warnings.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "session-info.txt"))
stopifnot(length(new_failures) == 0L, isTRUE(lab_result))
