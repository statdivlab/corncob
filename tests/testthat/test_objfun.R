library(corncob)
context("Test objfun failure")

set.seed(1)
seq_depth <- rpois(20, lambda = 10000)
my_counts <- rbinom(20, size = seq_depth, prob = 0.001) * 10
my_covariate <- cbind(rep(c(0,1), each = 10))
colnames(my_covariate) <- c("X1")

test_data <- data.frame("W" = my_counts, "M" = seq_depth, my_covariate)

should_break <- objfun(theta = c(1,1,1,-1000), W = my_counts, M = seq_depth,
                 X = cbind(1, my_covariate), X_star = cbind(1, my_covariate),
                 np = 2, npstar = 2,
                 link = "logit",
                 phi.link = "logit")

should_break_val <- objfun(theta = c(1,1,1,1000), W = my_counts, M = seq_depth,
                       X = cbind(1, my_covariate), X_star = cbind(1, my_covariate),
                       np = 2, npstar = 2,
                       link = "logit",
                       phi.link = "logit")


test_that("objfun breaks properly", {
  expect_equal(should_break$value, Inf)
  expect_equal(should_break_val$value, Inf)
})
