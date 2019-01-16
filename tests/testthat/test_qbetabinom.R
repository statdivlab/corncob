library(corncob)
context("Test qbetabinom")


test_that("qbetabinom works", {
  expect_is(qbetabinom(p = .5, size = 1000, mu = .1, phi = .1), "numeric")
  expect_is(qbetabinom(p = .001, size = 1000, mu = .001, phi = .1), "numeric")
})


test_that("HPDbetabinom works", {
  expect_is(HPDbetabinom(0, 10, .5, .1), "list")
})
