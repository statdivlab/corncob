library(corncob)
context("Test raotest")

set.seed(1)
seq_depth <- rpois(20, lambda = 10000)
my_counts <- rbinom(20, size = seq_depth, prob = 0.1)
my_covariate <- cbind(rep(c(0,1), each = 10), rep(c(0,1), 10))
colnames(my_covariate) <- c("X1", "X2")

test_data <- data.frame("W" = my_counts, "M" = seq_depth, my_covariate)

out1 <- bbdml(formula = cbind(W, M - W) ~ 1,
              phi.formula = ~ 1,
              data = test_data,
              link = "logit",
              phi.link = "logit",
              nstart = 1)

out2 <- bbdml(formula = cbind(W, M - W) ~ X1,
              phi.formula = ~ 1,
              data = test_data,
              link = "logit",
              phi.link = "logit",
              nstart = 1)

raotest(mod = out2, mod_null = out1)
lrtest(mod = out2, mod_null = out1)
waldt(out2)

test_that("raotest works", {
  expect_is(raotest(mod = out2, mod_null = out1), "numeric")
  expect_true(raotest(mod = out2, mod_null = out1) <= 1)
  expect_true(raotest(mod = out2, mod_null = out1) >= 0)
})

### confirm p-values are uniform under the null
### and high power under the alternative
# load_all()
# set.seed(3)
# nsims <- 100
# nn <- 20
# ps <- vector("numeric", length=nsims)
# for (i in 1:nsims) {
#
#   seq_depth <- rpois(nn, lambda = 10000)
#   my_covariate <- cbind(rep(c(0,1), each = nn/2))
#   colnames(my_covariate) <- c("X1")
#   probs <- rbeta(nn, 5, 5 + 0*c(my_covariate))
#   my_counts <- rbinom(nn, size = seq_depth, prob = probs)
#   test_data <- data.frame("W" = my_counts, "M" = seq_depth, my_covariate)
#
#   model_alt <- bbdml(formula = cbind(W, M - W) ~ X1,
#                phi.formula = ~ 1,
#                data = test_data,
#                link = "logit",
#                phi.link = "logit",
#                nstart = 1)
#
#   model_null <- bbdml(formula = cbind(W, M - W) ~ 1,
#                phi.formula = ~ 1,
#                data = test_data,
#                link = "logit",
#                phi.link = "logit",
#                nstart = 1)
#
#   ps[i] <- raotest(mod=model_alt, mod_null=model_null)
#
# }
# ps
# hist(ps)
# plot(nns, res[, 1])
# plot(nns, sqrt(nns)*res[, 1])
# plot(res[, 1], res_model[, 1], log="xy"); abline(0,1)
