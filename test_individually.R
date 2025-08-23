source('golden_tests/validators/R/validate_golden_tests.R')
suppressPackageStartupMessages(library(neuroim2))

test_files <- list.files("golden_tests/specs", pattern = "\\.xml$", 
                        recursive = TRUE, full.names = TRUE)

for (test_file in test_files) {
  cat("\nTesting:", basename(test_file), "... ")
  result <- tryCatch({
    validate_golden_test(test_file, verbose = FALSE)
  }, error = function(e) {
    list(success = FALSE, error = e$message)
  })
  
  if (result$success) {
    cat("✅ PASSED\n")
  } else {
    cat("❌ FAILED\n")
    if (!is.null(result$error)) {
      cat("  Error:", result$error, "\n")
    }
  }
}