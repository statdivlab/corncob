library(corncob)
context("Test links")

test_that("Link functions match identities", {
  expect_equal(logit(0.5), 0)
  expect_equal(invlogit(0), 0.5)
  expect_equal(fishZ(0), 0)
  expect_equal(invfishZ(0), 0)
  expect_equal(coth(1), 1/tanh(1))
})
