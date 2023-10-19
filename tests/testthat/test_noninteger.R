library(corncob)
context("Test noninteger Ws work")

set.seed(2)
seq_depth <- rpois(20, lambda = 10000) + rnorm(20)
my_counts <- rbinom(20, size = round(seq_depth), prob = 0.001) * 10 + rnorm(20)
my_covariate <- cbind(rep(c(0,1), each = 10))
colnames(my_covariate) <- c("X1")

test_data <- data.frame("W" = my_counts, "M" = seq_depth, my_covariate)

out_integer <- bbdml(formula = round(cbind(W, M - W)) ~ X1,
                     phi.formula = ~ X1,
                     data = test_data,
                     link = "logit",
                     phi.link = "logit",
                     nstart = 1)

test_that("Noninteger W's throw errors if not forced", {

  expect_error({
    out_noninteger_error <- bbdml(formula = cbind(W, M - W) ~ X1,
                            phi.formula = ~ X1,
                            data = test_data,
                            link = "logit",
                            phi.link = "logit",
                            nstart = 1)
  })
})

out_noninteger <- bbdml(formula = cbind(W, M - W) ~ X1,
                        phi.formula = ~ X1,
                        data = test_data,
                        link = "logit",
                        phi.link = "logit",
                        nstart = 1,
                        allow_noninteger=TRUE)

test_that("Noninteger W's work", {

  expect_is(summary(out_noninteger), "summary.bbdml")

  # Noninteger W's are close to integer W's
  expect_equal(out_noninteger$param, out_integer$param, tolerance=0.05)

  # Sandwich standard errors work for noninteger W's

  expect_is(sandSE(out_integer), "matrix")
})

