options(repos = c(CRAN = "https://cloud.r-project.org"))
stopifnot(rmarkdown::pandoc_available())
rcmdcheck::rcmdcheck(
  args = "--no-manual",
  build_args = "--no-manual",
  check_dir = "check",
  error_on = "warning"
)
