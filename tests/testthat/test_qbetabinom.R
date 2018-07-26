library(corncob)
context("Test qbetabinom")


test_that("qbetabinom works", {
  expect_is(qbetabinom(p = .5, size = 1000, mu = .1, phi = .1), "numeric")
  expect_is(qbetabinom(p = .001, size = 1000, mu = .001, phi = .1), "numeric")
})
