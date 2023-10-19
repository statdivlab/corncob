library(corncob)
context("Test that separation leads to expected results")

seq_depth <- rpois(20, lambda = 10000) + rnorm(20)
my_covariate <- cbind(rep(c(0,1), each = 10))
colnames(my_covariate) <- c("X1")
my_counts <- c(rbinom(10, size = round(seq_depth), prob = 0.001),
               rbinom(10, size = round(seq_depth), prob = 0.1))
test_data_sep <- data.frame("W" = my_counts, "M" = round(seq_depth), my_covariate)


test_that("Perfect separation but nonzero totals across one group", {
  out <- bbdml(formula = cbind(W, M - W) ~ X1,
               phi.formula = ~ X1,
               data = test_data_sep,
               link = "logit",
               phi.link = "logit",
               nstart = 1)

  expect_is(summary(out), "summary.bbdml")
})

my_counts <- c(rbinom(10, size = round(seq_depth), prob = 0),
               rbinom(10, size = round(seq_depth), prob = 0.1))
test_data_sep <- data.frame("W" = my_counts, "M" = round(seq_depth), my_covariate)

test_that("Perfect separation with zero totals across one group", {
  expect_warning(bbdml(formula = cbind(W, M - W) ~ X1,
                       phi.formula = ~ X1,
                       data = test_data_sep,
                       link = "logit",
                       phi.link = "logit",
                       nstart = 1))
})


my_counts <- c(rbinom(10, size = round(seq_depth), prob = 0),
               rbinom(10, size = round(seq_depth), prob = 0.1) + rnorm(10))
test_data_sep <- data.frame("W" = my_counts, "M" = round(seq_depth), my_covariate)

test_that("Perfect separation with zero totals across one group & continuous Ws", {
  expect_warning(bbdml(formula = cbind(W, M - W) ~ X1,
                       phi.formula = ~ X1,
                       data = test_data_sep,
                       link = "logit",
                       phi.link = "logit",
                       nstart = 1,
                       allow_noninteger=TRUE))
})
