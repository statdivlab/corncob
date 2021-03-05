#' Identify differentially-abundant and differentially-variable taxa
#'
#' @param formula an object of class \code{formula} without the response: a symbolic description of the model to be fitted to the abundance
#' @param phi.formula an object of class \code{formula} without the response: a symbolic description of the model to be fitted to the dispersion
#' @param formula_null Formula for mean under null, without response
#' @param phi.formula_null Formula for overdispersion under null, without response
#' @param data a data frame containing the OTU table, or \code{phyloseq} object containing the variables in the models
#' @param link link function for abundance covariates, defaults to \code{"logit"}
#' @param phi.link link function for dispersion covariates, defaults to \code{"logit"}
#' @param test Character. Hypothesis testing procedure to use. One of \code{"Wald"} or \code{"LRT"} (likelihood ratio test).
#' @param boot Boolean. Defaults to \code{FALSE}. Indicator of whether or not to use parametric bootstrap algorithm. (See \code{\link{pbWald}} and \code{\link{pbLRT}}).
#' @param B Optional integer. Number of bootstrap iterations. Ignored if \code{boot} is \code{FALSE}. Otherwise, defaults to \code{1000}.
#' @param sample_data Data frame or matrix. Defaults to \code{NULL}. If \code{data} is a data frame or matrix, this must be included as covariates/sample data.
#' @param taxa_are_rows Boolean. Optional. If \code{data} is a data frame or matrix, this indicates whether taxa are rows. Defaults to \code{TRUE}.
#' @param filter_discriminant Boolean. Defaults to \code{TRUE}. If \code{FALSE}, discriminant taxa will not be filtered out.
#' @param fdr_cutoff Integer. Defaults to \code{0.05}. Desired type 1 error rate
#' @param fdr Character. Defaults to \code{"fdr"}. False discovery rate control method, see \code{\link{p.adjust}} for more options.
#' @param full_output Boolean. Opetional. Defaults to \code{FALSE}. Indicator of whether to include full \code{bbdml} model output for all taxa.
#' @param inits Optional initializations for model fit using \code{formula} and \code{phi.formula} as rows of a matrix. Defaults to \code{NULL}.
#' @param inits_null Optional initializations for model fit using \code{formula_null} and \code{phi.formula_null} as rows of a matrix. Defaults to \code{NULL}.
#' @param try_only Optional numeric. Will try only the \code{try_only} taxa, specified either via numeric input or character taxa names. Useful for speed when troubleshooting. Defaults to \code{NULL}, testing all taxa.
#' @param ... Optional additional arguments for \code{\link{bbdml}}
#'
#' @details See package vignette for details and example usage. Make sure the number of columns in all of the initializations are correct! \code{inits} probably shouldn't match \code{inits_null}. To use a contrast matrix, see \code{\link{contrastsTest}}.
#'
#' @return An object of class \code{differentialTest}. List with elements \code{p} containing the p-values, \code{p_fdr} containing the p-values after false discovery rate control,  \code{significant_taxa} containing the taxa names of the statistically significant taxa, \code{significant_models} containing a list of the model fits for the significant taxa, \code{all_models} containing a list of the model fits for all taxa, \code{restrictions_DA} containing a list of covariates that were tested for differential abundance, \code{restrictions_DV} containing a list of covariates that were tested for differential variability, \code{discriminant_taxa_DA} containing the taxa for which at least one covariate associated with the abundance was perfectly discriminant, \code{discriminant_taxa_DV} containing the taxa for which at least one covariate associated with the dispersion was perfectly discriminant, \code{data} containing the data used to fit the models. If \code{full_output = TRUE}, it will also include \code{full_output}, a list of all model output from \code{bbdml}.
#'
#'
#' @examples
#' \donttest{
#' # phyloseq example
#' data(soil_phylo)
#' soil <- soil_phylo %>%
#' phyloseq::subset_samples(DayAmdmt %in% c(11,21)) %>%
#' phyloseq::tax_glom("Phylum")
#' da_analysis <- differentialTest(formula = ~ DayAmdmt,
#'                                 phi.formula = ~ DayAmdmt,
#'                                 formula_null = ~ 1,
#'                                 phi.formula_null = ~ DayAmdmt,
#'                                 test = "Wald", boot = FALSE,
#'                                 data = soil,
#'                                 fdr_cutoff = 0.05)
#' }
#' @export
differentialTest <- function(formula, phi.formula,
                             formula_null, phi.formula_null,
                             data,
                             link = "logit",
                             phi.link = "logit",
                             test,
                             boot = FALSE,
                             B = 1000,
                             sample_data = NULL,
                             taxa_are_rows = TRUE,
                             filter_discriminant = TRUE,
                             fdr_cutoff = 0.05,
                             fdr = "fdr",
                             full_output = FALSE,
                             inits = NULL,
                             inits_null = NULL,
                             try_only = NULL,
                             ...) {

  # Record call
  call <- match.call(expand.dots = TRUE)
  # Record mu link
  link <- match.arg(link, choices = "logit")
  # Record phi link
  phi.link <- match.arg(phi.link, choices = c("fishZ", "logit"))

  # Convert phyloseq objects
  if ("phyloseq" %in% class(data)) {
    # Set up response
    taxanames <- phyloseq::taxa_names(data)
  } else if (is.matrix(data) || is.data.frame(data)) {

    # use phyloseq
    OTU <- phyloseq::otu_table(data, taxa_are_rows = taxa_are_rows)

    # Make sample data
    sampledata <- phyloseq::sample_data(data.frame(
      sample_data,
      row.names = phyloseq::sample_names(OTU)
    ))

    # Make phyloseq object
    data <- phyloseq::phyloseq(OTU, sampledata)
    # Set up response
    taxanames <- phyloseq::taxa_names(data)
  } else {
    stop("Input must be either data frame, matrix, or phyloseq object!")
  }

  # Set up output
  pvals <- perfDisc_DA <- perfDisc_DV <- rep(NA, length(taxanames))
  model_summaries <- rep(list(NA), length(taxanames))
  full_outputs <- rep(list(NA), length(taxanames))
  # check to make sure inits is of the same length
  if (!is.null(inits)) {
    ncol1 <- ncol(stats::model.matrix(object = formula, data = data.frame(phyloseq::sample_data(data))))
    ncol2 <- ncol(stats::model.matrix(object = phi.formula, data = data.frame(phyloseq::sample_data(data))))
    if (length(inits) != ncol1 + ncol2) {
      stop("inits must match number of regression parameters in formula and phi.formula!")
    }
  }
  # inits_null_mu
  if (!is.null(inits_null)) {
    ncol1 <- ncol(stats::model.matrix(object = formula_null, data = data.frame(phyloseq::sample_data(data))))
    ncol2 <- ncol(stats::model.matrix(object = phi.formula_null, data = data.frame(phyloseq::sample_data(data))))
    if (length(inits_null) != ncol1 + ncol2) {
      stop("init_null must match number of regression parameters in formula_null and phi.formula_null!")
    }
  }

    restrict_ind <- 0

    if (is.null(try_only)) {
      try_only <- 1:length(taxanames)
    }

    if (is.character(try_only)) {
      try_only <- which(taxanames %in% try_only)
    }
    # Loop through OTU/taxa
    for (i in try_only) {

      # Subset data to only select that taxa
      data_i <- convert_phylo(data, select = taxanames[i])

      if (sum(data_i$W) == 0) {
        perfDisc_DA[i] <- TRUE
        perfDisc_DV[i] <- TRUE
      } else {
        # Update formula to match
        formula_i <- stats::update(formula, cbind(W, M - W) ~ .)
        formula_null_i <- stats::update(formula_null, cbind(W, M - W) ~ .)

        # Fit unrestricted model
        mod <- suppressWarnings(try(bbdml(formula = formula_i, phi.formula = phi.formula,
                                          data = data_i, link = link, phi.link = phi.link,
                                          inits = inits, ...), silent = TRUE))

        # Fit restricted model
        mod_null <- suppressWarnings(try(bbdml(formula = formula_null_i, phi.formula = phi.formula_null,
                                               data = data_i, link = link, phi.link = phi.link,
                                               inits = inits_null, ...), silent = TRUE))

        if (!("try-error" %in% c(class(mod), class(mod_null)))) {
          if (restrict_ind == 0) {
            restricts <- getRestrictionTerms(mod = mod, mod_null = mod_null)
            restrict_ind <- 1
          }
          # If both models fit, otherwise keep as NA
          model_summaries[[i]] <- suppressWarnings(summary(mod))
          if (full_output) {
            full_outputs[[i]] <- suppressWarnings(mod)
          }
          if (test == "Wald") {
            if (boot) {
              tmp <- try(pbWald(mod = mod, mod_null = mod_null, B = B), silent = TRUE)
              if (class(tmp) != "try-error") {
                pvals[i] <- tmp
              }
            } else {
              tmp <- try(waldchisq(mod = mod, mod_null = mod_null), silent = TRUE)
              if (class(tmp) != "try-error") {
                pvals[i] <- tmp
              }
            }
          } else if (test == "LRT") {
            if (boot) {
              tmp <- try(pbLRT(mod = mod, mod_null = mod_null, B = B), silent = TRUE)
              if (class(tmp) != "try-error") {
                pvals[i] <- tmp
              }
            } else {
              tmp <- try(lrtest(mod = mod, mod_null = mod_null), silent = TRUE)
              if (class(tmp) != "try-error") {
                pvals[i] <- tmp
              }
            }
          }
          perfDisc_DA[i] <- mod$sep_da
          perfDisc_DV[i] <- mod$sep_dv
        }
      }
    }

    ind_disc_da <- which(perfDisc_DA == TRUE)
    ind_disc_dv <- which(perfDisc_DV == TRUE)
    disc_vec_da <- taxanames[ind_disc_da]
    disc_vec_dv <- taxanames[ind_disc_dv]

    ind_disc <- union(ind_disc_da, ind_disc_dv)

    if (filter_discriminant && length(ind_disc) > 0) {
      # Want to keep same length, rest will ignore NAs
      pvals[ind_disc] <- NA
    }

    if (all(is.na(pvals))) {
      stop("All models failed to converge! \n
           If you are seeing this, it is likely that your model is overspecified. This occurs when your sample size is not large enough to estimate all the parameters of your model. This is most commonly due to categorical variables that include many categories. \n
           Alternatively, double-check your values for the arguments `link`, `phi.link`, and `method` to makes sure that they follow the specified options. \n
           To confirm you have fixed the issue, try running a model for a single taxon with bbdml.")
    }
    post_fdr <- stats::p.adjust(pvals, method = fdr)
    names(pvals) <- names(post_fdr) <- taxanames
    # Record significant taxa
    signif_vec <- taxanames[which(post_fdr < fdr_cutoff)]
    signif_models <- model_summaries[which(post_fdr < fdr_cutoff)]


    # restricts_mu <- setdiff(attr(terms(formula), "term.labels"),
    #                         attr(terms(formula_null), "term.labels"))
    # restricts_phi <- setdiff(attr(terms(phi.formula), "term.labels"),
    #                          attr(terms(phi.formula_null), "term.labels"))

    restricts_mu <- setdiff(colnames(stats::model.matrix(object = formula, data = data.frame(phyloseq::sample_data(data)))),
                            colnames(stats::model.matrix(object = formula_null, data = data.frame(phyloseq::sample_data(data)))))
    restricts_phi <- setdiff(colnames(stats::model.matrix(object = phi.formula, data = data.frame(phyloseq::sample_data(data)))),
                             colnames(stats::model.matrix(object = phi.formula_null, data = data.frame(phyloseq::sample_data(data)))))

    attr(restricts_mu, "index") <- restricts$mu
    attr(restricts_phi, "index") <- restricts$phi

    if (!full_output) {
      full_outputs <- NULL
    }


  structure(
    list("p" = pvals, "p_fdr" = post_fdr,
         "significant_taxa" = signif_vec,
         "significant_models" = signif_models,
         "all_models" =  model_summaries,
         "restrictions_DA" = restricts_mu,
         "restrictions_DV" = restricts_phi,
         "discriminant_taxa_DA" = disc_vec_da,
         "discriminant_taxa_DV" = disc_vec_dv,
         "data" = data,
         "full_output" = full_outputs),
    class = "differentialTest"
  )
}
