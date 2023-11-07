library(corncob)
suppressWarnings(library(phyloseq))
context("Test contrastsTest")

set.seed(1)
data(soil_phylo)
soil <- phyloseq::subset_samples(soil_phylo, DayAmdmt %in% c(11,21))
subsoil <- phyloseq::prune_taxa(x = soil, taxa = rownames(otu_table(soil))[301:325])
# including 7027 for all zeros
subsoil_disc <- phyloseq::prune_taxa(x = soil, taxa = rownames(otu_table(soil))[c(3001:3024, 7027)])

data(soil_phylum_small)
temp <- contrastsTest(formula = ~ DayAmdmt,
                             phi.formula = ~ DayAmdmt,
                             contrasts_DA = list("DayAmdmt21 - DayAmdmt11",
                                                 "DayAmdmt22 - DayAmdmt21"),
                             data = soil_phylum_small,
                             fdr_cutoff = 0.05)


test_that("contrastTest works", {
  expect_is(temp, "contrastsTest")
})
