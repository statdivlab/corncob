library(corncob)
context("Test qbetabinom")


test_that("qbetabinom works", {
  expect_is(qbetabinom(p = .5, M = 1000, mu = .1, phi = .1), "numeric")
  expect_is(qbetabinom(p = .001, M = 1000, mu = .001, phi = .1), "numeric")
})


test_that("HPDbetabinom works", {
  expect_is(HDIbetabinom(0, 10, .5, .1), "list")
})
