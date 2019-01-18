context("Test getRestrictionTerms")

set.seed(1)
seq_depth <- rpois(20, lambda = 10000)
my_counts <- rbinom(20, size = seq_depth, prob = 0.001) * 10
my_covariate <- cbind(rep(c(0,1), each = 10), rep(c(0,1), 10))
colnames(my_covariate) <- c("X1", "X2")

test_data <- data.frame("W" = my_counts, "M" = seq_depth, my_covariate)


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

out_interact <-  bbdml(formula = cbind(W, M - W) ~ X2*X1,
                       phi.formula = ~ X1,
                       data = test_data,
                       link = "logit",
                       phi.link = "logit",
                       nstart = 1)

out_noint <-  bbdml(formula = cbind(W, M - W) ~ X1,
                    phi.formula = ~ X1 - 1,
                    data = test_data,
                    link = "logit",
                    phi.link = "logit",
                    nstart = 1)

data(soil_phylo)
soil <- soil_phylo %>%
  phyloseq::subset_samples(DayAmdmt %in% c(11,21)) %>%
  phyloseq::tax_glom("Phylum")
mod1 <-  bbdml(formula = OTU.1 ~ Day*Plants,
               phi.formula = ~ Plants,
               data = soil)

mod2 <-  bbdml(formula = OTU.1 ~Day - 1,
               phi.formula = ~ Plants - 1,
               data = soil)

test_that("getRestrictionTerms works", {
  tmp <- corncob:::getRestrictionTerms(out,out_nullmu)
  expect_true(tmp$mu == 2)
  expect_null(tmp$phi)
  tmp <- corncob:::getRestrictionTerms(out,out_nullphi)
  expect_true(tmp$phi == 4)
  expect_null(tmp$mu)
  expect_true(attr(tmp$phi, "added"))
  expect_error(corncob:::getRestrictionTerms(out_nullmu, out_nullphi))
  expect_error(corncob:::getRestrictionTerms(out_nullphi, out_nullmu))
  expect_error(corncob:::getRestrictionTerms(out, restrictions = TRUE))
  tmp <- corncob:::getRestrictionTerms(out, restrictions.phi = 1)
  expect_true(tmp$phi == 3)
  expect_error(corncob:::getRestrictionTerms(out, restrictions.phi = TRUE))
  tmp <- corncob:::getRestrictionTerms(out_interact, out_noint)
  expect_equal(tmp$mu, c(2,4))
  expect_true(tmp$phi == 5)
  tmp <- corncob:::getRestrictionTerms(out, restrictions = 1)
  expect_equal(tmp$mu, 1)
  tmp <- corncob:::getRestrictionTerms(mod1,mod2)
  expect_equal(tmp$mu, c(1,3,4))
  expect_true(tmp$phi == 5)
  tmp <- corncob:::getRestrictionTerms(out, restrictions = "(Intercept)", restrictions.phi = "(Intercept)")
  expect_equal(tmp$mu, 1)
  expect_true(tmp$phi == 3)
})

