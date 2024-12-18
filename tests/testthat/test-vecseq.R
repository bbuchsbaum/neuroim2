library(purrr)
library(testthat)
library(assertthat)

gmask5 <- system.file("extdata", "global_mask_v4.nii", package="neuroim2")
mask <- system.file("extdata", "global_mask2.nii.gz", package="neuroim2")

context("neurovecseq")

test_that("can construct a NeuroVecSeq", {
  vec <- read_vec(gmask5)
  cvec <- concat(vec,vec,vec)

  vs <- NeuroVecSeq(vec,vec,vec)
  expect_equal(dim(cvec), dim(vs))
  expect_equal(space(cvec), space(vs))

})

test_that("can linearly index a NeuroVecSeq", {
  vec <- read_vec(gmask5)
  cvec <- concat(vec,vec,vec)

  vs <- NeuroVecSeq(vec,vec,vec)
  expect_equal(vs[1], cvec[1])
  for (i in 1:20) {
    idx <- sample(1:prod(dim(cvec)), 50)
    expect_equal(vs[idx], cvec[idx])
  }
})

test_that("can array index a NeuroVecSeq", {
  vec <- read_vec(gmask5)
  cvec <- concat(vec,vec,vec)

  vs <- NeuroVecSeq(vec,vec,vec)
  expect_equal(vs[1], cvec[1])

  expect_equal(vs[1,,,], cvec[1,,,])
  expect_equal(vs[1,2,,], cvec[1,2,,])
  expect_equal(vs[1,2,3,], cvec[1,2,3,])
  expect_equal(vs[1,2,3,4], cvec[1,2,3,4])
  expect_equal(vs[1:2,2,3,4], cvec[1:2,2,3,4])
  expect_equal(vs[1:2,2,3,3:4], cvec[1:2,2,3,3:4])
  expect_equal(vs[1:2,2,,3:4], cvec[1:2,2,,3:4])
  expect_equal(vs[1:2,,,3:4], cvec[1:2,,,3:4])
  expect_equal(vs[,,,3:4], cvec[,,,3:4])

})

test_that("can extract vectors from NeuroVecSeq", {
  vec <- read_vec(gmask5)
  cvec <- concat(vec,vec,vec)

  vs <- NeuroVecSeq(vec,vec,vec)
  v1 <- vs %>% vectors()
  v2 <- cvec %>% vectors()
  checkind <- sample(1:length(v1), 100)
  for (i in checkind) {
    expect_equal(v1[[i]], v2[[i]])
  }

})

test_that("can extract vols from NeuroVecSeq", {
  vec <- read_vec(gmask5)
  cvec <- concat(vec,vec,vec)

  vs <- NeuroVecSeq(vec,vec,vec)
  v1 <- vs %>% vols()
  v2 <- cvec %>% vols()
  for (i in length(v1)) {
    expect_equal(v1[[i]], v2[[i]])
  }

})


test_that("can map over vols from NeuroVecSeq", {
  vec <- read_vec(gmask5)
  cvec <- concat(vec,vec,vec)

  vs <- NeuroVecSeq(vec,vec,vec)
  out <- do.call(concat, vs %>% vols() %>% purrr::map( ~ .))
  #out <- do.call(concat, vlist)
  expect_true(all(as.vector(cvec) == as.vector(out)))
})

test_that("can subset a NeuroVecSeq", {
  vec <- read_vec(c(gmask5, gmask5, gmask5))
  vec5 <- sub_vector(vec, 1:4)
  vec2 <- read_vec(gmask5)
  expect_true(all(vec5[] == vec2[]))
})

test_that("can read a sparse NeuroVecSeq", {
  m <- read_vol(mask)
  vec <- read_vec(c(gmask5, gmask5, gmask5), mask=as.logical(m))
  vec5 <- sub_vector(vec, 1:4)
  vec2 <- read_vec(gmask5, mask=as.logical(m))

  ## TODO check that values are same
  expect_true(all(dim(vec5) == dim(vec2)))
})

test_that("can use a MappedNeuroVec as elements in a NeuroVecSeq", {
  vec <- read_vec(gmask5)
  vs <- NeuroVecSeq(vec,vec,vec)

  mvecs <- replicate(3, MappedNeuroVec(gmask5), simplify=FALSE)
  vs2 <- do.call(NeuroVecSeq, mvecs)


  v1 <- vs %>% vols()
  v2 <- vs2 %>% vols()

  for (i in 1:length(v1)) {
    expect_equal(v1[[i]], v2[[i]])
  }

  v1 <- vs %>% vectors()
  v2 <- vs2 %>% vectors()
  checkind <- sample(1:length(v1), 100)
  for (i in checkind) {
    expect_equal(v1[[i]], v2[[i]])
  }

})


