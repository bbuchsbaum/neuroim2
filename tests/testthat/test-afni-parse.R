library(testthat)
library(neuroim2)

# Expected values were produced by the earlier stringr-based implementation
# (str_split / str_trim / str_sub); the base-R helpers must match them exactly.

ns <- asNamespace("neuroim2")
parse_int <- get("parseIntAttribute", envir = ns)
parse_float <- get("parseFloatAttribute", envir = ns)
parse_string <- get("parseStringAttribute", envir = ns)
str_split1 <- get(".str_split1", envir = ns)
str_trim <- get(".str_trim", envir = ns)

test_that("integer and float attributes split on any whitespace", {
  expect_identical(parse_int("  1 2   3 "), c(1L, 2L, 3L))
  expect_identical(parse_int(c(" 64 64", "25  4\t")), c(64L, 64L, 25L, 4L))
  expect_identical(parse_int("7"), 7L)
  expect_identical(parse_int(""), NA_integer_)
  expect_identical(parse_int("\t-1  2 "), c(-1L, 2L))
  expect_identical(parse_float(" 1.5 -2e-3  4 "), c(1.5, -0.002, 4))
  expect_identical(parse_float(c("0.0 1", "2.25")), c(0, 1, 2.25))
  expect_identical(parse_float("3.5\r\n"), 3.5)
})

test_that("string attributes keep the stringr tilde semantics", {
  expect_identical(parse_string("'abc~def~"), c("abc", "def"))
  expect_identical(parse_string("'single~"), "single")
  expect_identical(parse_string("'no tilde"), "no tilde")
  expect_identical(parse_string("'a~~b~"), c("a", "", "b"))
  expect_identical(parse_string("''~"), "'")
  expect_identical(parse_string("'"), "")
  expect_identical(parse_string(c("'x~y~", "'z~")), c("x", "y"))
})

test_that("split and trim helpers match stringr edge cases", {
  expect_identical(str_split1("a~b~", "~"), c("a", "b", ""))
  expect_identical(str_split1("~a", "~"), c("", "a"))
  expect_identical(str_split1("", "~"), "")
  expect_identical(str_split1(NA_character_, "~"), NA_character_)
  expect_identical(str_trim(c(" a ", "\t b\n", " c ", NA, "")),
                   c("a", "b", "c", NA, ""))
})

test_that("read_afni_header parses a minimal header", {
  hdr <- tempfile(fileext = ".HEAD")
  on.exit(unlink(hdr))
  writeLines(c(
    "",
    "type = string-attribute",
    "name = BRICK_LABS",
    "count = 8",
    "'a~bb~",
    "",
    "type = integer-attribute",
    "name = DATASET_RANK",
    "count = 3",
    " 3 2 0",
    "",
    "type  =  float-attribute",
    "name = BRICK_FLOAT_FACS",
    "count = 2",
    " 0.5",
    " 2"
  ), hdr)
  h <- neuroim2:::read_afni_header(hdr)
  expect_named(h, c("BRICK_LABS", "DATASET_RANK", "BRICK_FLOAT_FACS"))
  expect_identical(h$BRICK_LABS$content, c("a", "bb"))
  expect_identical(h$BRICK_LABS$type, "string-attribute")
  expect_identical(h$DATASET_RANK$content, c(3L, 2L, 0L))
  expect_identical(h$DATASET_RANK$count, "3")
  expect_identical(h$BRICK_FLOAT_FACS$type, "float-attribute")
  expect_identical(h$BRICK_FLOAT_FACS$content, c(0.5, 2))
})
