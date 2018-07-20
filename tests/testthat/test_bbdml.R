library(corncob)
context("Test bbdml")

set.seed(1)
seq_depth <- rpois(20, lambda = 10000)
my_counts <- rbinom(20, size = seq_depth, prob = 0.001) * 10
my_covariate <- cbind(rep(c(0,1), each = 10))
colnames(my_covariate) <- c("X1")

test_data <- data.frame("W" = my_counts, "M" = seq_depth, my_covariate)

out_bfgs_with_inits <- bbdml(formula = cbind(W, M - W) ~ X1,
             phi.formula = ~ X1,
             data = test_data,
             link = "logit",
             phi.link = "logit",
             method = "BFGS",
             inits = rbind(c(1,1,1,1), c(2,2,2,2)))

test_that("bbdml with BFGS and inits works", {
  expect_is(out_bfgs_with_inits, "bbdml")
})
