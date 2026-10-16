# Run from the repository root in a fresh R process. Exercise each article's
# actual setup chunk with parameters that conflict with stale global state.
check_vignette_setup <- function() {
  old_wd <- setwd("vignettes")
  on.exit(setwd(old_wd), add = TRUE)
  knitr::opts_chunk$set(error = FALSE)

  assign("params", list(family = "red", preset = "homage"), envir = globalenv())
  article_params <- list(family = "lapis", preset = "study")
  expected_theme <- if (requireNamespace("albersdown", quietly = TRUE)) {
    albersdown::theme_albers(
      family = article_params$family, preset = article_params$preset
    )
  } else {
    neuroim2::theme_neuro()
  }

  files <- list.files(pattern = "[.]Rmd$")
  stopifnot(length(files) > 0L)
  for (file in files) {
    lines <- readLines(file, warn = FALSE)
    start <- grep("^```\\{r setup[,}]", lines)
    stopifnot(length(start) == 1L)
    end <- which(seq_along(lines) > start & lines == "```")[1L]
    stopifnot(!is.na(end))

    article_env <- new.env(parent = globalenv())
    article_env$params <- article_params
    knitr::knit(text = lines[start:end], envir = article_env, quiet = TRUE)

    if (!isTRUE(all.equal(ggplot2::theme_get(), expected_theme))) {
      stop(file, ": plotting theme did not use the article's parameters")
    }
    if (requireNamespace("albersdown", quietly = TRUE)) {
      stopifnot(identical(article_env$article_params, article_params))
    }
    stopifnot(
      exists("demo_mask", envir = article_env, inherits = FALSE),
      identical(environment(article_env$demo_mask), article_env),
      !exists("demo_mask", envir = globalenv(), inherits = FALSE),
      identical(globalenv()$params, list(family = "red", preset = "homage"))
    )
    cat("PASS", file, "isolated parameters and shared helpers\n")
  }
}

check_vignette_setup()
