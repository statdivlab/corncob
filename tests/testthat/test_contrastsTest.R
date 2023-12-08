library(corncob)
context("Test contrastsTest")

set.seed(1)

test_that("contrastTest works", {
  data(soil_phylum_small_contrasts)
  limma_install <- try(find.package("limma"), silent = TRUE)

  if (requireNamespace("phyloseq", quietly = TRUE)) {
    if (!(inherits(limma_install, "try-error"))) {
      temp <- contrastsTest(formula = ~ DayAmdmt,
                            phi.formula = ~ DayAmdmt,
                            contrasts_DA = list("DayAmdmt21 - DayAmdmt11",
                                                "DayAmdmt22 - DayAmdmt21"),
                            data = soil_phylum_small_contrasts,
                            fdr_cutoff = 0.05)
    }
    if (!(inherits(limma_install, "try-error"))) {
      expect_is(temp, "contrastsTest")
    }
  }
  else {
    expect_error(contrastsTest(formula = ~ DayAmdmt,
                               phi.formula = ~ DayAmdmt,
                               contrasts_DA = list("DayAmdmt21 - DayAmdmt11",
                                                   "DayAmdmt22 - DayAmdmt21"),
                               data = soil_phylum_small_contrasts,
                               fdr_cutoff = 0.05))
  }
})
