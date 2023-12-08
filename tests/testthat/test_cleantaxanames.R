library(corncob)
context("Test clean_taxa_names")

test_that("clean_taxa_names works", {
  if (requireNamespace("phyloseq", quietly = TRUE)) {

    data(soil_phylo)
    soil <- phyloseq::subset_samples(soil_phylo, DayAmdmt %in% c(11,21))

    tmp1 <- phyloseq::taxa_names(clean_taxa_names(soil))
    tmp2 <- phyloseq::taxa_names(clean_taxa_names(soil, name = "Seq"))

    expect_equal(length(tmp1), 7770)
    expect_equal(tmp1[1], "OTU1")
    expect_equal(tmp1[1234], "OTU1234")
    expect_equal(tmp2[2543], "Seq2543")
    expect_error(clean_taxa_names(c(1,2,3)))
  } else {
    expect_error(warn_phyloseq())
  }
})

