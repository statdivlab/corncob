library(corncob)
context("Test contrastsTest")

test_that("contrastTest works (non-phyloseq object)", {
  set.seed(1)
  limma_install <- try(find.package("limma"), silent = TRUE)
  data(soil_phylum_contrasts_sample)
  data(soil_phylum_contrasts_otu)

  if (!(inherits(limma_install, "try-error"))) {
    temp <- contrastsTest(formula = ~ DayAmdmt,
                          phi.formula = ~ DayAmdmt,
                          contrasts_DA = list("DayAmdmt21 - DayAmdmt11",
                                              "DayAmdmt22 - DayAmdmt21"),
                          data = soil_phylum_contrasts_otu,
                          sample_data = soil_phylum_contrasts_sample,
                          fdr_cutoff = 0.05)
    expect_is(temp, "contrastsTest")
  } else {
    expect_error(contrastsTest(formula = ~ DayAmdmt,
                               phi.formula = ~ DayAmdmt,
                               contrasts_DA = list("DayAmdmt21 - DayAmdmt11",
                                                   "DayAmdmt22 - DayAmdmt21"),
                               data = soil_phylum_contrasts_otu,
                               sample_data = soil_phylum_contrasts_sample,
                               fdr_cutoff = 0.05),
                 "If you would like to test contrasts, please install the `limma` package, available through Bioconductor.")
  }
})
