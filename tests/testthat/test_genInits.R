library(corncob)
context("Test genInits failure")


set.seed(1)
seq_depth <- rpois(20, lambda = 10000)
my_counts <- rbinom(20, size = seq_depth, prob = 0.001) * 10
my_covariate <- cbind(rep(c(0,1), each = 10))
colnames(my_covariate) <- c("X1")

test_data <- data.frame("W" = my_counts, "M" = seq_depth, my_covariate)

out <- genInits(nstart = 1, max.call = NULL, temperature = NULL,
         W = my_counts, M = seq_depth,
         X = cbind(1, my_covariate), X_star = cbind(1, my_covariate),
         np = 2, npstar = 2,
         link = "logit", phi.link = "logit", logpar = TRUE)

test_that("genInits works", {
  expect_is(out, "matrix")
})


test_that("genInits breaks properly", {
  expect_error(suppressWarnings(genInits(1, W = rep(c(1000,0), each = 10), M = seq_depth,
                        X = cbind(1, my_covariate), X_star = cbind(1, my_covariate),
                        np = 2, npstar = 2,
                        link = "logit", phi.link = "fishZ", logpar = TRUE, max.call = 100)))
})

