library(corncob)
library(phyloseq)
context("Test differentialTest")

set.seed(1)
data(soil_phylo)
soil <- phyloseq::subset_samples(soil_phylo, DayAmdmt %in% c(11,21))
subsoil <- prune_taxa(x = soil, taxa = rownames(otu_table(soil))[301:325])
temp <- differentialTest(formula = ~ Plants + DayAmdmt,
                         phi.formula = ~ Plants + DayAmdmt,
                         formula_null = ~ 1,
                         phi.formula_null = ~ 1,
                         data = subsoil,
                         inits = rbind(rep(.01, 6)),
                         inits_null_mu = rbind(rep(0.01, 4)),
                         inits_null_phi = rbind(rep(0.01, 4)))

temp_badinits1 <- differentialTest(formula = ~ Plants + DayAmdmt,
                         phi.formula = ~ Plants + DayAmdmt,
                         formula_null = ~ 1,
                         phi.formula_null = ~ 1,
                         data = subsoil,
                         inits = rbind(rep(Inf, 6)),
                         inits_null_mu = rbind(rep(0.01, 4)),
                         inits_null_phi = rbind(rep(0.01, 4)))

temp_badinits2 <- differentialTest(formula = ~ Plants + DayAmdmt,
                         phi.formula = ~ Plants + DayAmdmt,
                         formula_null = ~ 1,
                         phi.formula_null = ~ 1,
                         data = subsoil,
                         inits = rbind(rep(.01, 6)),
                         inits_null_mu = rbind(rep(Inf, 4)),
                         inits_null_phi = rbind(rep(Inf, 4)))

# Add this to cause some warnings and check those
temp_noinit_sing <- differentialTest(formula = ~ DayAmdmt,
                                phi.formula = ~ DayAmdmt,
                                formula_null = ~ 1,
                                phi.formula_null = ~ 1,
                                data = subsoil)

temp_noinit <- differentialTest(formula = ~ Plants + DayAmdmt,
                         phi.formula = ~ Plants + DayAmdmt,
                         formula_null = ~ 1,
                         phi.formula_null = ~ 1,
                         data = subsoil)

temp_sing <- differentialTest(formula = ~ DayAmdmt,
                              phi.formula = ~ DayAmdmt,
                              formula_null = ~ 1,
                              phi.formula_null = ~ 1,
                              data = subsoil,
                              inits = rbind(rep(.01, 4)))

temp_badinits3 <- differentialTest(formula = ~ DayAmdmt,
                              phi.formula = ~ DayAmdmt,
                              formula_null = ~ 1,
                              phi.formula_null = ~ 1,
                              data = subsoil,
                              inits = rbind(rep(Inf, 4)))

mydat <- phyloseq::get_taxa(subsoil)
mysampdat <- phyloseq::get_variable(subsoil)

temp_nonphylo <- differentialTest(formula = ~ DayAmdmt,
                                              phi.formula = ~ DayAmdmt,
                                              formula_null = ~ 1,
                                              phi.formula_null = ~ 1,
                                              data = mydat,
                                              sample_data = mysampdat,
                                              inits = rbind(rep(.01, 4)))


test_that("differentialTest works", {
  expect_is(temp, "list")
  expect_is(temp_sing, "list")
  expect_is(temp_nonphylo, "list")
  expect_is(temp_noinit, "list")
  expect_is(temp_noinit_sing, "list")
  expect_is(temp_badinits1, "list")
  expect_is(temp_badinits2, "list")
  expect_is(temp_badinits3, "list")
})

test_that("differentialTest works without phyloseq", {
  expect_true(all.equal(temp_sing, temp_nonphylo))
})


test_that("otu_to_taxonomy works", {
  expect_is(otu_to_taxonomy(temp$significant_taxa, soil_phylo), "character")
})

test_that("requires data frame, matrix, or phyloseq", {
  expect_error(differentialTest(formula = ~ DayAmdmt,
                                phi.formula = ~ DayAmdmt,
                                formula_null = ~ 1,
                                phi.formula_null = ~ 1,
                                data = c(1,2,3),
                                inits = rbind(rep(.01, 4))))
})

test_that("inits require correct length", {
  expect_error(differentialTest(formula = ~ DayAmdmt,
                                phi.formula = ~ DayAmdmt,
                                formula_null = ~ 1,
                                phi.formula_null = ~ 1,
                                data = subsoil,
                                inits = rbind(rep(.01, 6))))
  expect_error(differentialTest(formula = ~ DayAmdmt,
                                phi.formula = ~ DayAmdmt,
                                formula_null = ~ 1,
                                phi.formula_null = ~ 1,
                                data = subsoil,
                                inits_null_mu = rbind(rep(.01, 4))))
  expect_error(differentialTest(formula = ~ DayAmdmt,
                                phi.formula = ~ DayAmdmt,
                                formula_null = ~ 1,
                                phi.formula_null = ~ 1,
                                data = subsoil,
                                inits_null_phi = rbind(rep(.01, 4))))
})
