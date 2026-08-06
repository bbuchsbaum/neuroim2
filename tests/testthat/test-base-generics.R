library(testthat)
library(neuroim2)

# neuroim2 used to export an S4 generic named `scale` that had no methods at
# all, so attaching the package broke `base::scale()` for every input (GitHub
# issue #16). The generic is gone; these pin that base behaviour survives.
test_that("attaching neuroim2 does not mask base::scale()", {
  x <- c(1, 2, 3)
  expect_equal(as.numeric(scale(x)), as.numeric(base::scale(x)))

  m <- matrix(1:6, nrow = 3)
  expect_equal(scale(m), base::scale(m))
  expect_equal(scale(m, center = TRUE, scale = FALSE),
               base::scale(m, center = TRUE, scale = FALSE))

  # `scale` must resolve to the base function, not a neuroim2 generic.
  expect_false(isGeneric("scale", where = asNamespace("neuroim2")))
  expect_identical(environmentName(environment(scale)), "base")
})
