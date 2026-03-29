context("header API additional coverage")

library(neuroim2)

# ---- print.NeuroHeader optional branches ----

# Build a minimal NeuroHeader directly to exercise all print branches
make_header <- function(...) {
  base <- list(
    dim         = c(10L, 10L, 10L),
    pixdim      = c(0, 1, 1, 1, 2, 0, 0, 0),
    spacing     = c(1, 1, 1),
    origin      = c(0, 0, 0),
    trans       = diag(4),
    qform       = list(matrix = diag(4), code = NA_integer_),
    sform       = list(matrix = diag(4), code = NA_integer_),
    intent_code = NA_integer_,
    intent_name = NA_character_,
    descrip     = NA_character_,
    data_type   = NA_character_,
    bitpix      = NA_integer_,
    scl_slope   = 1,
    scl_inter   = 0,
    cal_min     = NA_real_,
    cal_max     = NA_real_,
    TR          = NA_real_,
    raw         = list()
  )
  extras <- list(...)
  for (nm in names(extras)) base[[nm]] <- extras[[nm]]
  class(base) <- "NeuroHeader"
  base
}

test_that("print.NeuroHeader shows sform code when present", {
  h <- make_header(sform = list(matrix = diag(4), code = 1L))
  out <- capture.output(print(h))
  expect_true(any(grepl("sform code", out)))
})

test_that("print.NeuroHeader shows qform code when present", {
  h <- make_header(qform = list(matrix = diag(4), code = 2L))
  out <- capture.output(print(h))
  expect_true(any(grepl("qform code", out)))
})

test_that("print.NeuroHeader shows TR when > 0", {
  h <- make_header(TR = 2.0, dim = c(10L, 10L, 10L, 20L))
  out <- capture.output(print(h))
  expect_true(any(grepl("TR", out)))
})

test_that("print.NeuroHeader shows intent_code when nonzero", {
  h <- make_header(intent_code = 2L)
  out <- capture.output(print(h))
  expect_true(any(grepl("Intent", out)))
})

test_that("print.NeuroHeader shows intent_name alongside intent_code", {
  h <- make_header(intent_code = 2L, intent_name = "TTEST")
  out <- capture.output(print(h))
  expect_true(any(grepl("TTEST", out)))
})

test_that("print.NeuroHeader shows descrip when present", {
  h <- make_header(descrip = "test description")
  out <- capture.output(print(h))
  expect_true(any(grepl("Description", out)))
})

test_that("print.NeuroHeader shows scaling when slope != 0 and != 1", {
  h <- make_header(scl_slope = 2.5, scl_inter = 10)
  out <- capture.output(print(h))
  expect_true(any(grepl("Scaling", out)))
})

test_that("print.NeuroHeader shows cal range when nonzero", {
  h <- make_header(cal_min = 0, cal_max = 100)
  out <- capture.output(print(h))
  expect_true(any(grepl("Cal range", out)))
})

test_that("print.NeuroHeader bitpix shown as 0 when NA", {
  h <- make_header(bitpix = NA_integer_, data_type = "FLOAT32")
  out <- capture.output(print(h))
  expect_true(any(grepl("0-bit", out)))
})

# ---- .build_neuro_header with minimal FileMetaInfo ----

test_that(".build_neuro_header from FileMetaInfo has all required fields", {
  f <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  skip_if(!file.exists(f), "extdata file not available")
  mi <- read_header(f)
  h <- neuroim2:::.build_neuro_header(mi)
  required <- c("dim", "pixdim", "spacing", "origin", "trans",
                "qform", "sform", "intent_code", "intent_name",
                "descrip", "data_type", "bitpix", "scl_slope",
                "scl_inter", "cal_min", "cal_max", "TR", "raw")
  expect_true(all(required %in% names(h)))
})

test_that("header returns NeuroHeader with expected class", {
  f <- system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  skip_if(!file.exists(f), "extdata file not available")
  h <- header(f)
  expect_true(inherits(h, "NeuroHeader"))
})
