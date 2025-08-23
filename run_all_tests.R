source('golden_tests/validators/R/validate_golden_tests.R')

# Suppress package loading messages
suppressPackageStartupMessages(library(neuroim2))

# Run all tests
results <- validate_all_golden_tests('golden_tests/specs', verbose = FALSE)

# Check for overall success
all_passed <- all(sapply(results, function(x) x$success))

if (all_passed) {
  cat("\n✅ ALL TESTS PASSED!\n")
} else {
  cat("\n❌ Some tests failed:\n")
  failed <- results[!sapply(results, function(x) x$success)]
  for (test in failed) {
    cat("  -", test$test_id, "\n")
  }
}