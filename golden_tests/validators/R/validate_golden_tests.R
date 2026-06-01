#' Golden Test Validator for R
#' 
#' This validator parses and executes golden test specifications in XML format,
#' validating numeric outputs against expected values with specified tolerances.

library(xml2)
library(testthat)

#' Parse a golden test XML file
#' 
#' @param xml_path Path to the golden test XML file
#' @return List containing test metadata, inputs, expected outputs, and R code
parse_golden_test <- function(xml_path) {
  doc <- read_xml(xml_path)
  
  # Extract metadata
  metadata <- list(
    id = xml_text(xml_find_first(doc, ".//d1:metadata/d1:id")),
    version = xml_text(xml_find_first(doc, ".//d1:metadata/d1:version")),
    description = xml_text(xml_find_first(doc, ".//d1:metadata/d1:description"))
  )
  
  # Extract semantic description
  semantic <- list(
    purpose = xml_text(xml_find_first(doc, ".//d1:semantic_description/d1:purpose")),
    algorithm = xml_text(xml_find_first(doc, ".//d1:semantic_description/d1:algorithm"))
  )
  
  # Extract R implementation
  r_code <- xml_text(xml_find_first(doc, ".//d1:implementations/d1:R"))
  
  # Extract numeric checks
  checks <- xml_find_all(doc, ".//d1:numeric_checks/d1:check")
  numeric_checks <- lapply(checks, function(check) {
    list(
      type = xml_text(xml_find_first(check, ".//d1:type")),
      name = xml_text(xml_find_first(check, ".//d1:name")),
      location = xml_text(xml_find_first(check, ".//d1:location")),
      expected = xml_text(xml_find_first(check, ".//d1:expected")),
      tolerance = as.numeric(xml_text(xml_find_first(check, ".//d1:tolerance"))),
      min = as.numeric(xml_text(xml_find_first(check, ".//d1:min"))),
      max = as.numeric(xml_text(xml_find_first(check, ".//d1:max")))
    )
  })
  
  list(
    metadata = metadata,
    semantic = semantic,
    r_code = r_code,
    numeric_checks = numeric_checks
  )
}

#' Perform a single numeric check
#' 
#' @param actual_value The actual value from test execution
#' @param check The check specification
#' @return TRUE if check passes, FALSE otherwise
perform_numeric_check <- function(actual_value, check) {
  expected_str <- trimws(check$expected)
  
  switch(check$type,
    "exact_value" = {
      expected <- as.numeric(strsplit(expected_str, " ")[[1]])
      tolerance <- ifelse(is.na(check$tolerance), 1e-10, check$tolerance)
      
      if (length(actual_value) != length(expected)) {
        return(FALSE)
      }
      all(abs(actual_value - expected) <= tolerance)
    },
    
    "dimension" = {
      expected <- as.numeric(strsplit(expected_str, " ")[[1]])
      identical(as.numeric(actual_value), expected)
    },
    
    "relative" = {
      expected <- as.numeric(expected_str)
      tolerance <- ifelse(is.na(check$tolerance), 0.01, check$tolerance)
      abs((actual_value - expected) / expected) <= tolerance
    },
    
    "range" = {
      all(actual_value >= check$min & actual_value <= check$max)
    },
    
    "statistical" = {
      # For statistical checks, implement specific logic based on the check name
      TRUE  # Placeholder
    },
    
    TRUE  # Default case
  )
}

#' Validate a golden test
#' 
#' @param xml_path Path to the golden test XML file
#' @param verbose Whether to print detailed output
#' @return List with test results
validate_golden_test <- function(xml_path, verbose = TRUE) {
  test_spec <- parse_golden_test(xml_path)
  
  if (verbose) {
    cat("\n=== Golden Test:", test_spec$metadata$id, "===\n")
    cat("Description:", test_spec$metadata$description, "\n")
  }
  
  # Create a new environment for test execution
  test_env <- new.env()
  
  # Execute the R code
  tryCatch({
    eval(parse(text = test_spec$r_code), envir = test_env)
  }, error = function(e) {
    return(list(
      test_id = test_spec$metadata$id,
      success = FALSE,
      error = paste("Code execution error:", e$message),
      checks = list()
    ))
  })
  
  # Perform numeric checks
  check_results <- lapply(test_spec$numeric_checks, function(check) {
    # Extract actual value from the test environment
    actual_expr <- parse(text = check$location)
    actual_value <- tryCatch({
      eval(actual_expr, envir = test_env)
    }, error = function(e) {
      if (verbose) {
        cat("  ERROR evaluating", check$location, ":", e$message, "\n")
      }
      return(NA)
    })
    
    if (all(is.na(actual_value))) {
      result <- list(
        name = check$name,
        passed = FALSE,
        message = paste("Could not evaluate:", check$location)
      )
    } else {
      passed <- tryCatch({
        perform_numeric_check(actual_value, check)
      }, error = function(e) {
        if (verbose) {
          cat("  ERROR in check:", e$message, "\n")
        }
        FALSE
      })
      
      result <- list(
        name = check$name,
        passed = passed,
        actual = actual_value,
        expected = check$expected,
        message = ifelse(passed, "PASS", 
                        paste("FAIL - Expected:", check$expected, 
                              "Got:", paste(actual_value, collapse = " ")))
      )
    }
    
    if (verbose) {
      cat(sprintf("  %-30s %s\n", check$name, result$message))
    }
    
    result
  })
  
  # Overall test result
  all_passed <- all(sapply(check_results, function(x) x$passed))
  
  if (verbose) {
    cat("\nOverall result:", ifelse(all_passed, "PASSED", "FAILED"), "\n")
  }
  
  list(
    test_id = test_spec$metadata$id,
    success = all_passed,
    checks = check_results
  )
}

#' Validate all golden tests in a directory
#' 
#' @param test_dir Directory containing golden test XML files
#' @param pattern File pattern to match (default: "\\.xml$")
#' @param verbose Whether to print detailed output
#' @return Summary of test results
validate_all_golden_tests <- function(test_dir, pattern = "\\.xml$", verbose = TRUE) {
  # Find all XML files recursively
  test_files <- list.files(test_dir, pattern = pattern, 
                          recursive = TRUE, full.names = TRUE)
  
  if (length(test_files) == 0) {
    cat("No golden test files found in", test_dir, "\n")
    return(NULL)
  }
  
  cat("Found", length(test_files), "golden test files\n")
  
  # Run all tests
  results <- lapply(test_files, function(test_file) {
    validate_golden_test(test_file, verbose = verbose)
  })
  
  # Summary
  passed <- sum(sapply(results, function(x) x$success))
  failed <- length(results) - passed
  
  cat("\n=== SUMMARY ===\n")
  cat("Total tests:", length(results), "\n")
  cat("Passed:", passed, "\n")
  cat("Failed:", failed, "\n")
  
  # Return detailed results
  invisible(results)
}

#' Run golden tests as part of testthat
#' 
#' @param test_dir Directory containing golden test XML files
#' @export
test_golden <- function(test_dir = "golden_tests/specs") {
  test_that("all golden tests pass", {
    results <- validate_all_golden_tests(test_dir, verbose = FALSE)
    
    for (result in results) {
      for (check in result$checks) {
        expect_true(check$passed, 
                   label = paste(result$test_id, "-", check$name),
                   info = check$message)
      }
    }
  })
}

# Example usage:
if (FALSE) {
  # Validate a single test
  validate_golden_test("golden_tests/specs/core/spatial_reference/neurospace_construction.xml")
  
  # Validate all tests
  validate_all_golden_tests("golden_tests/specs")
  
  # Use with testthat
  test_golden()
}