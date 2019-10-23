library(corncob)
context("Test waldtest")

set.seed(1)
seq_depth <- rpois(20, lambda = 10000)
my_counts <- rbinom(20, size = seq_depth, prob = 0.001) * 10
my_covariate <- cbind(rep(c(0,1), each = 10))
colnames(my_covariate) <- c("X1")

test_data <- data.frame("W" = my_counts, "M" = seq_depth, my_covariate)

my_counts_bad <- my_counts
my_counts_bad[1:19] <- 0
my_covariate <- cbind(rep(c(0,1), each = 10), rep(c(1,0), each = 10))
colnames(my_covariate) <- c("X1", "X2")
test_data_bad <- data.frame("W" = my_counts_bad, "M" = seq_depth, my_covariate)

out <- bbdml(formula = cbind(W, M - W) ~ X1,
             phi.formula = ~ X1,
             data = test_data,
             link = "logit",
             phi.link = "logit",
             nstart = 1)
out_nullmu <- bbdml(formula = cbind(W, M - W) ~ 1,
                    phi.formula = ~ X1,
                    data = test_data,
                    link = "logit",
                    phi.link = "logit",
                    nstart = 1)
out_nullphi <- bbdml(formula = cbind(W, M - W) ~ X1,
                    phi.formula = ~ 1,
                    data = test_data,
                    link = "logit",
                    phi.link = "logit",
                    nstart = 1)

out_nullboth <- bbdml(formula = cbind(W, M - W) ~ 1,
                      phi.formula = ~ 1,
                      data = test_data,
                      link = "logit",
                      phi.link = "logit",
                      nstart = 1)

out_bad <- bbdml(formula = cbind(W, M - W) ~ X1+X2,
                 phi.formula = ~ X1+X2,
                 data = test_data_bad,
                 link = "logit",
                 phi.link = "fishZ",
                 nstart = 1)

test_that("waldt works", {
  expect_is(waldt(out), "matrix")
})

test_that("waldchisq works", {
  expect_is(waldchisq(out, out_nullmu), "numeric")
  expect_is(waldchisq(out, out_nullphi), "numeric")
  tmp1 <- waldchisq(mod = out, restrictions = 2, restrictions.phi = 2)
  tmp2 <- waldchisq(mod = out, mod_null = out_nullboth)
  tmp3 <- waldchisq(mod = out, restrictions = "X1", restrictions.phi = "X1")
  tmp4 <- waldchisq(mod = out, restrictions = "X1", restrictions.phi = 2)
  tmp5 <- waldchisq(mod = out, restrictions = 2, restrictions.phi = "X1")
  expect_is(tmp1, "numeric")
  expect_equal(tmp1, tmp2)
  expect_equal(tmp1, tmp3)
  expect_equal(tmp1, tmp4)
  expect_equal(tmp1, tmp5)
  expect_true(is.na(waldchisq(out, restrictions = 5)))
})

test_that("waldtest can break", {
  expect_error(waldt(c(1,2,3)))
  expect_warning(waldt(out_bad))
})

test_that("waldchisq_test works", {
  expect_error(corncob:::waldchisq_test(out, restrictions = integer(0)))
  expect_error(corncob:::waldchisq_test(out, restrictions = TRUE))
  expect_error(corncob:::waldchisq_test(out, restrictions.phi = TRUE))
  expect_is(corncob:::waldchisq_test(out, restrictions = "X1", restrictions.phi = "X1"), "numeric")
  expect_is(corncob:::waldchisq_test(out, restrictions.phi = "X1"), "numeric")
  expect_is(corncob:::waldchisq_test(out, restrictions = "X1"), "numeric")
  expect_error(corncob:::waldchisq_test(out_bad, restrictions = 2))
})
