library(corncob)
suppressWarnings(library(phyloseq))
context("Test bbdml")

set.seed(1)
seq_depth <- rpois(20, lambda = 10000)
my_counts <- rbinom(20, size = seq_depth, prob = 0.001) * 10
my_covariate <- cbind(rep(c(0,1), each = 10))
colnames(my_covariate) <- c("X1")

test_data <- data.frame("W" = my_counts, "M" = seq_depth, my_covariate)


my_counts_bad <- my_counts
my_counts_bad[1:19] <- 0
my_covariate <- cbind(rep(c(0,1), each = 10), c(rep(0,5), rep(1,15)), c(rep(0,4), rep(1,16)))
colnames(my_covariate) <- c("X1", "X2", "X3")
test_data_bad <- data.frame("W" = my_counts_bad, "M" = seq_depth, my_covariate)

test_small <- test_data[c(1,11),]

test_that("overspecified model fails", {
  expect_error(bbdml(formula = cbind(W, M - W) ~ X1,
                     phi.formula = ~ X1,
                     data = test_small,
                     link = "logit",
                     phi.link = "logit",
                     method = "trust",
                     inits = rbind(c(1,1,1,1), c(2,2,2,2))))
})

out_bfgs_inits_num <- bbdml(formula = cbind(W, M - W) ~ X1,
             phi.formula = ~ X1,
             data = test_data,
             link = "logit",
             phi.link = "logit",
             method = "BFGS",
             numerical = TRUE,
             inits = rbind(c(1,1,1,1), c(2,2,2,2)))

out_trust <- bbdml(formula = cbind(W, M - W) ~ X1,
                   phi.formula = ~ X1,
                   data = test_data,
                   link = "logit",
                   phi.link = "logit",
                   method = "trust",
                   inits = rbind(c(1,1,1,1), c(2,2,2,2)))

out_bad_init <- suppressWarnings(bbdml(formula = cbind(W, M - W) ~ X1,
                   phi.formula = ~ X1,
                   data = test_data,
                   link = "logit",
                   phi.link = "logit",
                   method = "trust",
                   inits = cbind(1,-1000,1,1000)))
test_that("bbdml with BFGS, inits, and numerical works", {
  expect_is(out_bfgs_inits_num, "bbdml")
  expect_is(out_trust, "bbdml")
  expect_is(out_bad_init, "bbdml")
})

test_that("bad init gives warning", {
  expect_warning(bbdml(formula = cbind(W, M - W) ~ X1,
                       phi.formula = ~ X1,
                       data = test_data,
                       link = "logit",
                       phi.link = "logit",
                       method = "trust",
                       inits = cbind(1,-1000,1,1000)))
})

test_that("checking for perfectly discriminant", {
  expect_warning(tmp <- bbdml(formula = cbind(W, M - W) ~ X1,
                       phi.formula = ~ X1,
                       data = test_data_bad,
                       link = "logit",
                       phi.link = "logit",
                       nstart = 1))
  expect_output(expect_warning(print(tmp)))
  expect_output(expect_warning(print(summary(tmp))))
  expect_true(tmp$sep_da) # confirm separation detected
  expect_true(tmp$sep_dv)
  expect_false(out_bad_init$sep_da) # sanity check comparison, no separation

  expect_warning(bbdml(formula = cbind(W, M - W) ~ 1,
                       phi.formula = ~ X1,
                       data = test_data_bad,
                       link = "logit",
                       phi.link = "logit",
                       nstart = 1))
  expect_warning(bbdml(formula = cbind(W, M - W) ~ X1,
                       phi.formula = ~ 1,
                       data = test_data_bad,
                       link = "logit",
                       phi.link = "logit",
                       nstart = 1))
  expect_warning(bbdml(formula = cbind(W, M - W) ~ X1-1,
                       phi.formula = ~ 1,
                       data = test_data_bad,
                       link = "logit",
                       phi.link = "logit",
                       nstart = 1))
  expect_warning(bbdml(formula = cbind(W, M - W) ~ 1,
                       phi.formula = ~ X1-1,
                       data = test_data_bad,
                       link = "logit",
                       phi.link = "logit",
                       nstart = 1))
})

data(soil_phylo)
soil <- phyloseq::subset_samples(soil_phylo, DayAmdmt %in% c(11,21))

out_phylo <- bbdml(formula = OTU.4 ~ 1,
             phi.formula = ~ 1,
             data = soil)


test_that("bbdml works with phyloseq object", {
  expect_is(out_phylo, "bbdml")
})

