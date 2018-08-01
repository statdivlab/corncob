library(corncob)
library(phyloseq)
context("Test differentialTest")

set.seed(1)
data(soil_phylo)
soil <- phyloseq::subset_samples(soil_phylo, DayAmdmt %in% c(11,21))

subsoil <- prune_taxa(x = soil, taxa = rownames(otu_table(soil))[301:325])
temp <- differentialTest(formula = ~ DayAmdmt,
                         phi.formula = ~ DayAmdmt,
                         formula_null = ~ 1,
                         phi.formula_null = ~ 1,
                         data = subsoil)

test_that("differentialTest works", {
  expect_is(temp, "list")
})

test_that("differentialTest breaks without phyloseq", {
  expect_error(suppressWarnings(differentialTest(formula = ~ DayAmdmt,
                                phi.formula = ~ DayAmdmt,
                                formula_null = ~ 1,
                                phi.formula_null = ~ 1,
                                data = c(123))))
})

test_that("differentialTest breaks with too many covariates", {
  expect_error(suppressWarnings(differentialTest(formula = ~ DayAmdmt + Plants,
                                                 phi.formula = ~ DayAmdmt,
                                                 formula_null = ~ 1,
                                                 phi.formula_null = ~ 1,
                                                 data = subsoil)))
})

test_that("otu_to_taxonomy works", {
  expect_is(otu_to_taxonomy(temp$DA, soil_phylo), "character")
})
