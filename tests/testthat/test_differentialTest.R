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
                         data = subsoil, boot = FALSE, test = "LRT",
                         inits = rbind(rep(.01, 6)),
                         inits_null = rbind(rep(0.01, 2)))

temp2 <- differentialTest(formula = ~ Plants + DayAmdmt,
                         phi.formula = ~ Plants + DayAmdmt,
                         formula_null = ~ 1,
                         phi.formula_null = ~ 1,
                         data = subsoil, boot = FALSE, test = "LRT",
                         filter_discriminant = FALSE,
                         inits = rbind(rep(.01, 6)),
                         inits_null = rbind(rep(0.01, 2)))

temp_badinits1 <- differentialTest(formula = ~ Plants + DayAmdmt,
                         phi.formula = ~ Plants + DayAmdmt,
                         formula_null = ~ 1,
                         phi.formula_null = ~ 1,
                         data = subsoil, boot = FALSE, test = "LRT",
                         inits = rbind(rep(Inf, 6)),
                         inits_null = rbind(rep(0.01, 2)))

temp_badinits2 <- differentialTest(formula = ~ Plants + DayAmdmt,
                         phi.formula = ~ Plants + DayAmdmt,
                         formula_null = ~ 1,
                         phi.formula_null = ~ 1,
                         data = subsoil, boot = FALSE, test = "LRT",
                         inits = rbind(rep(.01, 6)),
                         inits_null = rbind(rep(Inf, 2)))

# Add this to cause some warnings and check those
temp_noinit_sing <- differentialTest(formula = ~ DayAmdmt,
                                phi.formula = ~ DayAmdmt,
                                formula_null = ~ 1,
                                phi.formula_null = ~ 1, boot = FALSE, test = "LRT",
                                data = subsoil)

temp_noinit <- differentialTest(formula = ~ Plants + DayAmdmt,
                         phi.formula = ~ Plants + DayAmdmt,
                         formula_null = ~ 1,
                         phi.formula_null = ~ 1, boot = FALSE, test = "LRT",
                         data = subsoil)

temp_sing <- differentialTest(formula = ~ DayAmdmt,
                              phi.formula = ~ DayAmdmt,
                              formula_null = ~ 1,
                              phi.formula_null = ~ 1,
                              data = subsoil, boot = FALSE, test = "LRT",
                              inits = rbind(rep(.01, 4)))

temp_badinits3 <- differentialTest(formula = ~ DayAmdmt,
                              phi.formula = ~ DayAmdmt,
                              formula_null = ~ 1,
                              phi.formula_null = ~ 1,
                              data = subsoil, boot = FALSE, test = "LRT",
                              inits = rbind(rep(Inf, 4)))

mydat <- phyloseq::get_taxa(subsoil)
mysampdat <- phyloseq::get_variable(subsoil)

temp_nonphylo <- differentialTest(formula = ~ DayAmdmt,
                                              phi.formula = ~ DayAmdmt,
                                              formula_null = ~ 1,
                                              phi.formula_null = ~ 1,
                                              data = mydat, boot = FALSE, test = "LRT",
                                              sample_data = mysampdat,
                                              inits = rbind(rep(.01, 4)))

temp_wald <- differentialTest(formula = ~ Plants + DayAmdmt,
                         phi.formula = ~ Plants + DayAmdmt,
                         formula_null = ~ 1,
                         phi.formula_null = ~ 1,
                         data = subsoil, boot = FALSE, test = "Wald",
                         inits = rbind(rep(.01, 6)),
                         inits_null = rbind(rep(0.01, 2)))

temp_pblrt <- differentialTest(formula = ~ Plants + DayAmdmt,
                         phi.formula = ~ Plants + DayAmdmt,
                         formula_null = ~ 1,
                         phi.formula_null = ~ 1,
                         data = subsoil, boot = TRUE, B = 5, test = "LRT",
                         inits = rbind(rep(.01, 6)),
                         inits_null = rbind(rep(0.01, 2)))

temp_pbwald <- differentialTest(formula = ~ Plants + DayAmdmt,
                               phi.formula = ~ Plants + DayAmdmt,
                               formula_null = ~ 1,
                               phi.formula_null = ~ 1,
                               data = subsoil, boot = TRUE, B = 5, test = "Wald",
                               inits = rbind(rep(.01, 6)),
                               inits_null = rbind(rep(0.01, 2)))


test_that("differentialTest works", {
  expect_is(temp, "differentialTest")
  expect_is(temp2, "differentialTest")
  expect_is(temp_wald, "differentialTest")
  expect_is(temp_pbwald, "differentialTest")
  expect_is(temp_pblrt, "differentialTest")
  expect_is(temp, "differentialTest")
  expect_is(temp_sing, "differentialTest")
  expect_is(temp_nonphylo, "differentialTest")
  expect_is(temp_noinit, "differentialTest")
  expect_is(temp_noinit_sing, "differentialTest")
  expect_is(temp_badinits1, "differentialTest")
  expect_is(temp_badinits2, "differentialTest")
  expect_is(temp_badinits3, "differentialTest")
})


test_that("differentialTest S3 methods", {
  expect_is(plot(temp), "ggplot")
  expect_is(plot(temp, level = c("Order", "Class")), "ggplot")
  expect_is(plot(temp, level = "Kingdom"), "ggplot")
  expect_null(print(temp))
})

test_that("differentialTest works without phyloseq", {
  expect_true(all.equal(temp_sing$p, temp_nonphylo$p))
})


test_that("otu_to_taxonomy works", {
  expect_is(otu_to_taxonomy(temp$significant_taxa, soil_phylo), "character")
})

test_that("requires data frame, matrix, or phyloseq", {
  expect_error(differentialTest(formula = ~ DayAmdmt,
                                phi.formula = ~ DayAmdmt,
                                formula_null = ~ 1,
                                phi.formula_null = ~ 1,
                                data = c(1,2,3), boot = FALSE, test = "LRT",
                                inits = rbind(rep(.01, 4))))
})

test_that("inits require correct length", {
  expect_error(differentialTest(formula = ~ DayAmdmt,
                                phi.formula = ~ DayAmdmt,
                                formula_null = ~ 1,
                                phi.formula_null = ~ 1,
                                data = subsoil, boot = FALSE, test = "LRT",
                                inits = rbind(rep(.01, 6))))
  expect_error(differentialTest(formula = ~ DayAmdmt,
                                phi.formula = ~ DayAmdmt,
                                formula_null = ~ 1,
                                phi.formula_null = ~ 1,
                                data = subsoil, boot = FALSE, test = "LRT",
                                inits_null = rbind(rep(.01, 4))))
})
