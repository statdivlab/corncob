library(corncob)
context("Test differentialTest")

data("soil_phylum_small_sample")
data("soil_phylum_small_otu")

temp <- differentialTest(formula = ~ Plants + DayAmdmt,
                         phi.formula = ~ Plants + DayAmdmt,
                         formula_null = ~ 1,
                         phi.formula_null = ~ 1,
                         data = soil_phylum_small_otu,
                         sample_data = soil_phylum_small_sample,
                         boot = FALSE, test = "LRT",
                         inits = rbind(rep(.01, 6)),
                         inits_null = rbind(rep(0.01, 2)))

test_that("differentialTest works for data frames", {
  expect_is(temp, "differentialTest")
})





