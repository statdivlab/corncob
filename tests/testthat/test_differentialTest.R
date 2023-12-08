library(corncob)
context("Test differentialTest")

set.seed(1)
if (requireNamespace("phyloseq", quietly = TRUE)) {
  data(soil_phylo)
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
}





