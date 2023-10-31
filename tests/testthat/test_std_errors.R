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
  expect_is(sand_vcov(out), "matrix")
})


test_that("Sandwich SEs are close to model-based SEs", {

  ## n is large here, and model is correct, so SEs should align (here, to within 5%)

  sandwich_ses <- sand_vcov(out) %>% diag %>% sqrt %>% unname
  model_ses <- coef(summary(out))[,2] %>% unname
  relative_differences <- abs(sandwich_ses - model_ses)/model_ses
  relative_differences

  expect_true(all(relative_differences < 0.05))

})


### confirm they scale correctly with 1/sqrt(n) and match model-based SEs -- yes, they do!


# set.seed(2)
# nsims <- 50
# nns <- 2*round(seq(from = 30, to = 500, length.out = nsims))
# res <- matrix(NA, ncol = 4, nrow = nsims)
# res_model <- matrix(NA, ncol = 4, nrow = nsims)
# for (i in 1:nsims) {
#   nn <- nns[i]
#   seq_depth <- rpois(nn, lambda = 10000)
#   probs <- rbeta(nn, 1, 1)
#   my_counts <- rbinom(nn, size = seq_depth, prob = probs)
#   my_covariate <- cbind(rep(c(0,1), each = nn/2))
#   colnames(my_covariate) <- c("X1")
#
#   test_data <- data.frame("W" = my_counts, "M" = seq_depth, my_covariate)
#
#   out <- bbdml(formula = cbind(W, M - W) ~ X1,
#                phi.formula = ~ X1,
#                data = test_data,
#                link = "logit",
#                phi.link = "logit",
#                nstart = 1)
#
#
#   res[i, ] <- sand_vcov(out) %>% diag %>% sqrt %>% unname
#   res_model[i, ] <- coef(summary(out))[,2] %>% unname
#
# }
# plot(nns, res[, 1])
# plot(nns, res_model[, 1])
# plot(nns, sqrt(nns)*res[, 1])
# plot(nns, sqrt(nns)*res_model[, 1])
# plot(res[, 1], res_model[, 1], log="xy"); abline(0,1)

