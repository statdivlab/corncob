library(corncob)
context("Test raotest")

test_that("raotest controls T1ER and is correlated w other p-values", {

  set.seed(5)
  nsims <- 50
  nn <- 100
  ps_rao <- vector("numeric", length=nsims)
  ps_rao_robust <- vector("numeric", length=nsims)
  ps_lrt <- vector("numeric", length=nsims)
  for (i in 1:nsims) {
    seq_depth <- rpois(nn, lambda = 10000)
    my_covariate <- cbind(rep(c(0,1), each = nn/2))
    my_covariate2 <- cbind(runif(nn))
    colnames(my_covariate) <- c("X1")
    colnames(my_covariate2) <- c("X2")
    probs <- rbeta(nn, 5, 5 + 0*c(my_covariate) + 1*c(my_covariate2))
    my_counts <- rbinom(nn, size = seq_depth, prob = probs)
    test_data <- data.frame("W" = my_counts, "M" = seq_depth, "X1" = my_covariate, "X2" = my_covariate2)

    model_alt <- bbdml(formula = cbind(W, M - W) ~ X1 + X2,
                       phi.formula = ~ X1,
                       data = test_data,
                       link = "logit",
                       phi.link = "logit",
                       nstart = 5)

    model_null <- bbdml(formula = cbind(W, M - W) ~ X2,
                        phi.formula = ~ 1,
                        data = test_data,
                        link = "logit",
                        phi.link = "logit",
                        nstart = 5)

    ps_rao[i] <- raotest(mod=model_alt, mod_null=model_null)
    model_alt$robust <- TRUE
    model_null$robust <- TRUE
    ps_rao_robust[i] <- raotest(mod=model_alt, mod_null=model_null)
    ps_lrt[i] <- lrtest(mod=model_alt, mod_null=model_null)

  }

  expect_true(all(ps_lrt <= 1))
  expect_true(all(ps_lrt >= 0))
  expect_true(all(ps_rao <= 1))
  expect_true(all(ps_rao >= 0))
  expect_true(all(ps_rao_robust <= 1))
  expect_true(all(ps_rao_robust >= 0))

  ### error rate control
  expect_true(mean(ps_rao < 0.05) <= 0.05 + 1.96*sqrt(0.05*0.95/nsims))
  expect_true(mean(ps_rao_robust < 0.05) <= 0.05 + 1.96*sqrt(0.05*0.95/nsims))
  expect_true(mean(ps_lrt < 0.05) <= 0.05 + 1.96*sqrt(0.05*0.95/nsims))

  ## correlations
  expect_true(cor(ps_rao, ps_lrt) > 0.9)
  expect_true(cor(ps_rao, ps_rao_robust) > 0.95)

})


test_that("raotest and lrtest is decently powered", {

  set.seed(6)
  nsims <- 50
  nn <- 200
  ps_rao <- vector("numeric", length=nsims)
  ps_rao_robust <- vector("numeric", length=nsims)
  ps_lrt <- vector("numeric", length=nsims)
  for (i in 1:nsims) {
    seq_depth <- rpois(nn, lambda = 10000)
    my_covariate <- cbind(rep(c(0,1), each = nn/2))
    my_covariate2 <- cbind(runif(nn))
    colnames(my_covariate) <- c("X1")
    colnames(my_covariate2) <- c("X2")
    probs <- rbeta(nn, 5, 5 + 0.5*c(my_covariate) + 1*c(my_covariate2))
    my_counts <- rbinom(nn, size = seq_depth, prob = probs)
    test_data <- data.frame("W" = my_counts, "M" = seq_depth, "X1" = my_covariate, "X2" = my_covariate2)

    model_alt <- bbdml(formula = cbind(W, M - W) ~ X1 + X2,
                       phi.formula = ~ 1,
                       data = test_data,
                       link = "logit",
                       phi.link = "logit",
                       nstart = 5)

    model_null <- bbdml(formula = cbind(W, M - W) ~ X2,
                        phi.formula = ~ 1,
                        data = test_data,
                        link = "logit",
                        phi.link = "logit",
                        nstart = 5)

    ps_rao[i] <- raotest(mod=model_alt, mod_null=model_null)
    model_alt$robust <- TRUE
    model_null$robust <- TRUE
    ps_rao_robust[i] <- raotest(mod=model_alt, mod_null=model_null)
    ps_lrt[i] <- lrtest(mod=model_alt, mod_null=model_null)

  }

  ### power. 10% is low, but double a "noninformative" test
  expect_true(mean(ps_rao < 0.05) > 0.1)
  expect_true(mean(ps_rao_robust < 0.05) > 0.1)
  expect_true(mean(ps_lrt < 0.05) > 0.1)

  ## correlation with LRT
  expect_true(cor(ps_rao, ps_lrt) > 0.9)
  expect_true(cor(ps_rao, ps_rao_robust) > 0.9)

})
