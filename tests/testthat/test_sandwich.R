library(corncob)
context("Test sandwich standard errors")

set.seed(1)
nn <- 2000
seq_depth <- rpois(nn, lambda = 10000)
probs <- rbeta(nn, 1, 1)
my_counts <- rbinom(nn, size = seq_depth, prob = probs)
my_covariate <- cbind(rep(c(0,1), each = nn/2))
colnames(my_covariate) <- c("X1")

test_data <- data.frame("W" = my_counts, "M" = seq_depth, my_covariate)

out <- bbdml(formula = cbind(W, M - W) ~ X1,
             phi.formula = ~ X1,
             data = test_data,
             link = "logit",
             phi.link = "logit",
             nstart = 1)

test_that("Sandwich standard errors work", {
  expect_is(sandSE(out), "matrix")
})


test_that("Sandwich SEs are close to model-based SEs", {

  ## n is large here, and model is correct, so SEs should align (here, to within 5%)

  sandwich_ses <- sandSE(out) %>% diag %>% sqrt %>% unname
  model_ses <- coef(summary(out))[,2] %>% unname
  relative_differences <- abs(sandwich_ses - model_ses)/model_ses
  relative_differences

  expect_true(all(relative_differences < 0.05))

})
