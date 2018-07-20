library(corncob)
context("Test Hessian and Score")

set.seed(1)
seq_depth <- rpois(20, lambda = 10000)
my_counts <- rbinom(20, size = seq_depth, prob = 0.001) * 10
my_covariate <- cbind(rep(c(0,1), each = 10))
colnames(my_covariate) <- c("X1")

test_data <- data.frame("W" = my_counts, "M" = seq_depth, my_covariate)

out <- bbdml(formula = cbind(W, M - W) ~ X1,
             phi.formula = ~ X1,
             data = test_data,
             link = "logit",
             phi.link = "logit",
             nstart = 1)

out_fish <- bbdml(formula = cbind(W, M - W) ~ X1,
                  phi.formula = ~ X1,
                  data = test_data,
                  link = "logit",
                  phi.link = "fishZ",
                  nstart = 1)

test_that("Analytic hessian works", {
  expect_is(hessian(out), "matrix")
  expect_is(hessian(out_fish), "matrix")
})

test_that("Numerical hessian works", {
  expect_is(hessian(out, numerical = TRUE), "matrix")
})

test_that("Analytic gradient works", {
  expect_is(score(out), "numeric")
  expect_is(score(out_fish), "numeric")
})

test_that("Numeric gradient works", {
  expect_is(score(out, numerical = TRUE), "numeric")
})

test_that("Subject-specific gradient works", {
  expect_is(score(out, forHess = TRUE), "matrix")
})
