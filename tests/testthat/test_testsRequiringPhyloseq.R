library(corncob)
context("Tests that use phyloseq")

if (requireNamespace("phyloseq", quietly = TRUE)) {
  # make phyloseq objects from data
  soil_phylo <- phyloseq::phyloseq(phyloseq::sample_data(soil_phylo_sample),
                                   phyloseq::otu_table(soil_phylo_otu, taxa_are_rows = TRUE),
                                   phyloseq::tax_table(soil_phylo_taxa))
  soil_phylum_small <- phyloseq::phyloseq(phyloseq::sample_data(soil_phylum_small_sample),
                                          phyloseq::otu_table(soil_phylum_small_otu, taxa_are_rows = TRUE))
  soil_phylum_small_contrasts <- phyloseq::phyloseq(phyloseq::sample_data(soil_phylum_contrasts_sample),
                                                    phyloseq::otu_table(soil_phylum_contrasts_otu, taxa_are_rows = TRUE))

  # tests for bbdml()
  test_that("bbdml works with phyloseq object", {
    out_phylo <- bbdml(formula = Proteobacteria ~ 1,
                       phi.formula = ~ 1,
                       data = soil_phylum_small)
    expect_is(out_phylo, "bbdml")
  })

  # tests for clean_taxa_names()
  test_that("clean_taxa_names works", {
    soil <- phyloseq::subset_samples(soil_phylo, DayAmdmt %in% c(11,21))

    tmp1 <- phyloseq::taxa_names(clean_taxa_names(soil))
    tmp2 <- phyloseq::taxa_names(clean_taxa_names(soil, name = "Seq"))

    expect_equal(length(tmp1), 7770)
    expect_equal(tmp1[1], "OTU1")
    expect_equal(tmp1[1234], "OTU1234")
    expect_equal(tmp2[2543], "Seq2543")
    expect_error(clean_taxa_names(c(1,2,3)))
  })

  # tests for contrastsTest()
  set.seed(1)
  limma_install <- try(find.package("limma"), silent = TRUE)
  if (!(inherits(limma_install, "try-error"))) {
    temp <- contrastsTest(formula = ~ DayAmdmt,
                          phi.formula = ~ DayAmdmt,
                          contrasts_DA = list("DayAmdmt21 - DayAmdmt11",
                                              "DayAmdmt22 - DayAmdmt21"),
                          data = soil_phylum_small_contrasts,
                          fdr_cutoff = 0.05)
  }

  test_that("contrastTest works", {
    if (!(inherits(limma_install, "try-error"))) {
      expect_is(temp, "contrastsTest")
    } else {
      expect_error(contrastsTest(formula = ~ DayAmdmt,
                                 phi.formula = ~ DayAmdmt,
                                 contrasts_DA = list("DayAmdmt21 - DayAmdmt11",
                                                     "DayAmdmt22 - DayAmdmt21"),
                                 data = soil_phylum_small_contrasts,
                                 fdr_cutoff = 0.05),
                   "If you would like to test contrasts, please install the `limma` package, available through Bioconductor.")
    }
  })

  # tests for getRestrictionTerms()
  test_that("getRestrictionTerms works for phyloseq object", {
    soil <- soil_phylo %>%
      phyloseq::subset_samples(DayAmdmt %in% c(11,21)) %>%
      phyloseq::tax_glom("Phylum")
    mod1 <-  bbdml(formula = OTU.1 ~ Day*Plants,
                   phi.formula = ~ Plants,
                   data = soil)

    mod2 <-  bbdml(formula = OTU.1 ~Day - 1,
                   phi.formula = ~ Plants - 1,
                   data = soil)

    tmp <- corncob:::getRestrictionTerms(mod1,mod2)
    expect_equal(tmp$mu, c(1,3,4))
    expect_true(tmp$phi == 5)
  })

  # tests for differentialTest()
  set.seed(1)
  soil <- phyloseq::subset_samples(soil_phylo, DayAmdmt %in% c(11,21))
  subsoil <- phyloseq::prune_taxa(x = soil, taxa = rownames(phyloseq::otu_table(soil))[301:325])
  # including 7027 for all zeros
  subsoil_disc <- phyloseq::prune_taxa(x = soil, taxa = rownames(phyloseq::otu_table(soil))[c(3001:3024, 7027)])

  subsoil_nonint <- subsoil
  phyloseq::otu_table(subsoil_nonint) <- phyloseq::otu_table(subsoil_nonint) + rexp(1)

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

  temp3 <- suppressWarnings(
    differentialTest(formula = ~ Plants + DayAmdmt,
                     phi.formula = ~ Plants + DayAmdmt,
                     formula_null = ~ 1,
                     phi.formula_null = ~ 1,
                     data = subsoil_disc, boot = FALSE, test = "LRT",
                     inits = rbind(rep(.01, 6)),
                     inits_null = rbind(rep(0.01, 2)))
    # Record of expected warnings below #

    # Warning in print.bbdml(mod) :
    # This model is based on a discriminant taxa.
    # You may see NAs in the model summary because Wald testing is invalid.
    # Likelihood ratio testing can be used, but valid standard errors cannot be calculated.

    # Warning in waldt(object) :
    # Singular Hessian! Cannot calculate p-values in this setting.
  )

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

  temp_pbwald_robust <- differentialTest(formula = ~ Plants + DayAmdmt,
                                         phi.formula = ~ Plants + DayAmdmt,
                                         formula_null = ~ 1,
                                         phi.formula_null = ~ 1,
                                         data = subsoil, boot = TRUE, B = 5, test = "Wald",
                                         inits = rbind(rep(.01, 6)),
                                         inits_null = rbind(rep(0.01, 2)),
                                         robust = TRUE)


  test_that("differentialTest works", {
    expect_is(temp, "differentialTest")
    expect_is(temp2, "differentialTest")
    expect_is(temp3, "differentialTest")
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
    expect_is(temp_pbwald_robust, "differentialTest")
  })


  test_that("differentialTest S3 methods", {
    expect_is(plot(temp), "ggplot")
    expect_is(plot(temp, level = c("Order", "Class")), "ggplot")
    expect_is(plot(temp, level = "Kingdom"), "ggplot")
    expect_output(expect_null(print(temp)))
  })

  test_that("differentialTest works without phyloseq", {
    expect_true(all.equal(temp_sing$p, temp_nonphylo$p))
  })


  test_that("otu_to_taxonomy works", {
    expect_is(otu_to_taxonomy(temp$significant_taxa, soil_phylo), "character")
    expect_error(otu_to_taxonomy(temp$significant_taxa, mysampdat))
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

  test_that("try_only works", {
    expect_is(differentialTest(formula = ~ DayAmdmt,
                               phi.formula = ~ DayAmdmt,
                               formula_null = ~ 1,
                               phi.formula_null = ~ 1,
                               data = subsoil, boot = FALSE, test = "LRT",
                               try_only = 1:2), "differentialTest")
  })

  test_that("overspecification error message", {
    expect_error(differentialTest(formula = ~ Plants*Day*Amdmt,
                                  phi.formula = ~ Plants*Day*Amdmt,
                                  formula_null = ~ 1,
                                  phi.formula_null = ~ 1,
                                  data = subsoil, boot = FALSE, test = "LRT",
                                  try_only = 1:2))
  })

  test_that("differentialTest does NAs correctly", {
    expect_equal(length(temp$p), 25)
    expect_equal(length(temp$p_fdr), 25)
    expect_equal(length(temp3$p), 25)
    expect_equal(length(temp3$p_fdr), 25)
  })



  test_that("differentialTest and non integers", {

    ## non integers should give error
    expect_error(
      temp5 <- differentialTest(formula = ~ Plants + DayAmdmt,
                                phi.formula = ~ Plants + DayAmdmt,
                                formula_null = ~ 1,
                                phi.formula_null = ~ 1,
                                data = subsoil_nonint,
                                boot = FALSE, test = "LRT")
    )

    # Warning with Wald, *non-robust*, boot = F
    expect_warning(
      temp6 <- differentialTest(formula = ~ Plants + DayAmdmt,
                                phi.formula = ~ Plants + DayAmdmt,
                                formula_null = ~ 1,
                                phi.formula_null = ~ 1,
                                data = subsoil_nonint,
                                allow_noninteger = TRUE,
                                boot = FALSE,
                                test = "Wald",
                                robust = FALSE)
    )

    # Everything goes through with Wald, robust, boot = F
    expect_silent(
      temp7 <- differentialTest(formula = ~ Plants + DayAmdmt,
                                phi.formula = ~ Plants + DayAmdmt,
                                formula_null = ~ 1,
                                phi.formula_null = ~ 1,
                                data = subsoil_nonint,
                                allow_noninteger = TRUE,
                                boot = FALSE,
                                test = "Wald",
                                robust = TRUE)
    )

    # Expect message with noninteger data with boot = T
    # because bootstrap is only parametric, and doesn't make sense with EE
    expect_message(
      # using invisible and capture.output to avoid printed output giving error from a failed
      # run of testing testing function
      invisible(utils::capture.output(temp8 <- differentialTest(formula = ~ Plants + DayAmdmt,
                                                                phi.formula = ~ Plants + DayAmdmt,
                                                                formula_null = ~ 1,
                                                                phi.formula_null = ~ 1,
                                                                data = subsoil_nonint,
                                                                boot = TRUE,
                                                                B = 5,
                                                                test = "Wald",
                                                                robust = TRUE,
                                                                allow_noninteger = TRUE)))
    )

    expect_silent(
      temp11 <- differentialTest(formula = ~ Plants + DayAmdmt,
                                 phi.formula = ~ Plants + DayAmdmt,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ 1,
                                 data = subsoil_nonint,
                                 allow_noninteger = TRUE,
                                 test = "Wald",
                                 robust=TRUE)
    )
    expect_warning(
      temp12 <- differentialTest(formula = ~ Plants + DayAmdmt,
                                 phi.formula = ~ Plants + DayAmdmt,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ 1,
                                 data = subsoil_nonint,
                                 allow_noninteger = TRUE,
                                 test = "Wald",
                                 robust=FALSE)
    )

    # Amy needs to think about whether this is allowed!!
    expect_error(
      temp13 <- differentialTest(formula = ~ Plants + DayAmdmt,
                                 phi.formula = ~ Plants + DayAmdmt,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ 1,
                                 data = subsoil_nonint,
                                 allow_noninteger = TRUE,
                                 test = "LRT",
                                 boot=FALSE)
    )

    expect_error(
      temp14 <- differentialTest(formula = ~ Plants + DayAmdmt,
                                 phi.formula = ~ Plants + DayAmdmt,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ 1,
                                 data = subsoil_nonint,
                                 allow_noninteger = TRUE,
                                 test = "LRT",
                                 B=3,
                                 boot=TRUE)
    )


    expect_silent(
      temp15 <- differentialTest(formula = ~ Plants + DayAmdmt,
                                 phi.formula = ~ Plants,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ 1,
                                 data = subsoil,
                                 allow_noninteger = TRUE,
                                 test = "Rao",
                                 robust=TRUE,
                                 B=3,
                                 boot=TRUE)
    )

    expect_silent(
      temp16 <- differentialTest(formula = ~ Plants + DayAmdmt,
                                 phi.formula = ~ Plants,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ 1,
                                 data = subsoil,
                                 allow_noninteger = TRUE,
                                 test = "Rao",
                                 robust=FALSE,
                                 B=3,
                                 boot=TRUE)
    )

    expect_silent(
      temp17 <- differentialTest(formula = ~ Plants + DayAmdmt,
                                 phi.formula = ~ Plants,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ 1,
                                 data = subsoil_nonint,
                                 allow_noninteger = TRUE,
                                 test = "Rao",
                                 robust=TRUE,
                                 B=3,
                                 boot=FALSE)
    )

    # expect warning when non-integers are present and robust = FALSE from Rao test
    expect_warning(
      temp18 <- differentialTest(formula = ~ Plants + DayAmdmt,
                                 phi.formula = ~ Plants,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ 1,
                                 data = subsoil_nonint,
                                 allow_noninteger = TRUE,
                                 test = "Rao",
                                 robust=FALSE,
                                 B=3,
                                 boot=FALSE)
    )

    expect_false(identical(temp17$p, temp18$p))
  })


} else {
  expect_error(warn_phyloseq())
}
