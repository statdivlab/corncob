#' Identify differentially-abundant and differentially-variable taxa using contrasts
#'
#' @param formula an object of class \code{formula} without the response: a symbolic description of the model to be fitted to the abundance
#' @param phi.formula an object of class \code{formula} without the response: a symbolic description of the model to be fitted to the dispersion
#' @param contrasts_DA List. Optional. Constructs a contrast matrix. List elements should be characters specifying contrasts in the parameters within \code{formula}. Note that this is only available with \code{"Wald"} value for \code{test}. Must include at least one of \code{contrasts_DA} or \code{contrasts_DV}.
#' @param contrasts_DV List. Optional. Constructs a contrast matrix. List elements should be characters specifying contrasts in the parameters within \code{phi.formula}. Note that this is only available with \code{"Wald"} value for \code{test}. Must include at least one of \code{contrasts_DA} or \code{contrasts_DV}.
#' @param data a data frame containing the OTU table, or \code{phyloseq} object containing the variables in the models
#' @param link link function for abundance covariates, defaults to \code{"logit"}
#' @param phi.link link function for dispersion covariates, defaults to \code{"logit"}
#' @param sample_data Data frame or matrix. Defaults to \code{NULL}. If \code{data} is a data frame or matrix, this must be included as covariates/sample data.
#' @param taxa_are_rows Boolean. Optional. If \code{data} is a data frame or matrix, this indicates whether taxa are rows. Defaults to \code{TRUE}.
#' @param filter_discriminant Boolean. Defaults to \code{TRUE}. If \code{FALSE}, discriminant taxa will not be filtered out.
#' @param fdr_cutoff Integer. Defaults to \code{0.05}. Desired type 1 error rate
#' @param fdr Character. Defaults to \code{"fdr"}. False discovery rate control method, see \code{\link{p.adjust}} for more options.
#' @param inits Optional initializations for model fit using \code{formula} and \code{phi.formula} as rows of a matrix. Defaults to \code{NULL}.
#' @param try_only Optional numeric. Will try only the \code{try_only} taxa, specified either via numeric input or character taxa names. Useful for speed when troubleshooting. Defaults to \code{NULL}, testing all taxa.
#' @param ... Optional additional arguments for \code{\link{bbdml}}
#'
#' @details This function uses contrast matrices to test for differential abundance and differential variability using a Wald-type chi-squared test. To use a formula implementation, see \code{\link{differentialTest}}.
#'
#' @return An object of class \code{contrastsTest}. List with elements \code{p} containing the p-values for each contrast, \code{p_fdr} containing the p-values after false discovery rate control,  \code{significant_taxa} containing the taxa names of the statistically significant taxa,  \code{contrasts_DA} containing the contrast matrix for parameters associated with the abundance, \code{contrasts_DV} containing the contrast matrix for parameters associated with the dispersion, \code{discriminant_taxa_DA} containing the taxa for which at least one covariate associated with the abundance was perfectly discriminant, \code{discriminant_taxa_DV} containing the taxa for which at least one covariate associated with the dispersion was perfectly discriminant, and \code{data} containing the data used to fit the models.

#' @examples
#'
#' # data frame example
#' # note that this function will only run if the `limma` package is installed
#' limma_install <- try(find.package("limma"), silent = TRUE)
#' if (!(inherits(limma_install, "try-error"))) {
#'   data(soil_phylum_contrasts_sample)
#'   data(soil_phylum_contrasts_otu)
#'   da_analysis <- contrastsTest(formula = ~ DayAmdmt,
#'                               phi.formula = ~ DayAmdmt,
#'                               contrasts_DA = list("DayAmdmt21 - DayAmdmt11",
#'                                                    "DayAmdmt22 - DayAmdmt21"),
#'                               data = soil_phylum_contrasts_otu,
#'                                sample_data = soil_phylum_contrasts_sample,
#'                                fdr_cutoff = 0.05,
#'                               try_only = 1:5)
#' }
#'
#' # phyloseq example (only run if you have phyloseq installed)
#' \dontrun{
#' contrasts_phylo <- phyloseq::phyloseq(phyloseq::sample_data(soil_phylum_contrasts_sample),
#' phyloseq::otu_table(soil_phylum_contrasts_otu, taxa_are_rows = TRUE))
#' da_analysis <- contrastsTest(formula = ~ DayAmdmt,
#'                              phi.formula = ~ DayAmdmt,
#'                              contrasts_DA = list("DayAmdmt21 - DayAmdmt11",
#'                                                  "DayAmdmt22 - DayAmdmt21"),
#'                              data = contrasts_phylo,
#'                              fdr_cutoff = 0.05,
#'                              try_only = 1:5)
#' }
#'
#' @export
contrastsTest <- function(formula, phi.formula,
                          contrasts_DA = NULL,
                          contrasts_DV = NULL,
                          data,
                          link = "logit",
                          phi.link = "logit",
                          sample_data = NULL,
                          taxa_are_rows = TRUE,
                          filter_discriminant = TRUE,
                          fdr_cutoff = 0.05,
                          fdr = "fdr",
                          inits = NULL,
                          try_only = NULL,
                          ...) {

  if (is.null(contrasts_DA) && is.null(contrasts_DV)) {
    stop("Must include at least one of contrasts_DA or contrasts_DV!")
  }

  limma_install <- try(find.package("limma"), silent = TRUE)
  if (inherits(limma_install, "try-error")) {
    {stop("If you would like to test contrasts, please install the `limma` package, available through Bioconductor.")}
  }

  num_DA <- length(contrasts_DA)
  num_DV <- length(contrasts_DV)

  # Record call
  call <- match.call(expand.dots = TRUE)
  # Record mu link
  link <- match.arg(link, choices = "logit")
  # Record phi link
  phi.link <- match.arg(phi.link, choices = c("fishZ", "logit"))

  # Convert phyloseq objects
  if ("phyloseq" %in% class(data)) {
    if (requireNamespace("phyloseq", quietly = TRUE)) {
      # Set up response
      taxanames <- phyloseq::taxa_names(data)
      sample_data <- data.frame(phyloseq::sample_data(data))
    } else {
      warn_phyloseq()
    }
  } else if (is.matrix(data) || is.data.frame(data)) {

    # # use phyloseq
    # OTU <- phyloseq::otu_table(data, taxa_are_rows = taxa_are_rows)
    #
    # # Make sample data
    # sampledata <- phyloseq::sample_data(data.frame(
    #   sample_data,
    #   row.names = phyloseq::sample_names(OTU)
    # ))
    #
    # # Make phyloseq object
    # data <- phyloseq::phyloseq(OTU, sampledata)
    # # Set up response
    # taxanames <- phyloseq::taxa_names(data)

    if (taxa_are_rows) {
      data <- t(data)
    }
    taxanames <- colnames(data)
    M <- rowSums(data)

  } else {
    stop("Input must be either data frame, matrix, or phyloseq object!")
  }

  # Set up output
  pvals <- perfDisc_DA <- perfDisc_DV <- matrix(NA, nrow = length(taxanames),
                                                ncol = num_DA + num_DV)
  #model_summaries <- rep(list(NA), length(taxanames))
  # check to make sure inits is of the same length
  if (!is.null(inits)) {
    ncol1 <- ncol(stats::model.matrix(object = formula, data = sample_data))
    ncol2 <- ncol(stats::model.matrix(object = phi.formula, data = sample_data))
    if (length(inits) != ncol1 + ncol2) {
      stop("inits must match number of regression parameters in formula and phi.formula!")
    }
  }

  if (is.null(try_only)) {
    try_only <- 1:length(taxanames)
  }

  if (is.character(try_only)) {
    try_only <- which(taxanames %in% try_only)
  }
  # Loop through OTU/taxa
  for (i in try_only) {

    # Subset data to only select that taxa
    if ("phyloseq" %in% class(data)) {
      data_i <- convert_phylo(data, select = taxanames[i])
    } else {
      response_i <- data.frame(W = data[, taxanames[i]], M = M)
      data_i <- cbind(response_i, sample_data)
    }

    if (sum(data_i$W) == 0) {
      perfDisc_DA[i] <- TRUE
      perfDisc_DV[i] <- TRUE
    } else {
      # Update formula to match
      formula_i <- stats::update(formula, cbind(W, M - W) ~ .)

      # Fit unrestricted model
      mod <- suppressWarnings(try(bbdml(formula = formula_i, phi.formula = phi.formula,
                                        data = data_i, link = link, phi.link = phi.link,
                                        inits = inits, ...), silent = TRUE))

      if (!inherits(mod, "try-error")) {
        # If both models fit, otherwise keep as NA
        #model_summaries[[i]] <- suppressWarnings(summary(mod))

        if (num_DA > 0) {
          for (contr in 1:num_DA) {
            tmp <- try(waldchisq(mod = mod, contrasts_DA = contrasts_DA[[contr]],
                                 contrasts_DV = NULL), silent = TRUE)
            if (!inherits(tmp, "try-error")) {
              pvals[i, contr] <- tmp
            }
            perfDisc_DA[i, contr] <- mod$sep_da
            perfDisc_DV[i, contr] <- mod$sep_dv
          }
        }

        if (num_DV > 0) {
          for (contr in 1:num_DV) {
            tmp <- try(waldchisq(mod = mod, contrasts_DA = NULL,
                                 contrasts_DV = contrasts_DV[[contr]]), silent = TRUE)
            if (!inherits(tmp, "try-error")) {
              pvals[i, contr + num_DA] <- tmp
            }
            perfDisc_DA[i, contr + num_DA] <- mod$sep_da
            perfDisc_DV[i, contr + num_DA] <- mod$sep_dv
          }
        }

      }
    }
  }

  disc_vec_da <- disc_vec_dv <- signif_vec <- rep(list(NA), num_DA + num_DV)
  post_fdr <- matrix(NA, nrow = length(taxanames), ncol = num_DA + num_DV)

  for (contr in 1:(num_DA + num_DV)) {
    ind_disc_da <- which(perfDisc_DA[, contr] == TRUE)
    ind_disc_dv <- which(perfDisc_DV[, contr] == TRUE)
    disc_vec_da[[contr]] <- taxanames[ind_disc_da]
    disc_vec_dv[[contr]] <- taxanames[ind_disc_dv]

    ind_disc <- union(ind_disc_da, ind_disc_dv)

    if (filter_discriminant && length(ind_disc) > 0) {
      # Want to keep same length, rest will ignore NAs
      pvals[ind_disc, contr] <- NA
    }

    if (all(is.na(pvals[, contr]))) {
      next
    }
    post_fdr[, contr] <- stats::p.adjust(pvals[, contr], method = fdr)
    names(pvals[, contr]) <- names(post_fdr[, contr]) <- taxanames
    # Record significant taxa
    signif_vec[[contr]] <- taxanames[which(post_fdr[, contr] < fdr_cutoff)]
    #signif_models <- model_summaries[which(post_fdr < fdr_cutoff)]
  }




  structure(
    list("p" = pvals, "p_fdr" = post_fdr,
         "significant_taxa" = signif_vec,
         #"significant_models" = signif_models,
         #"all_models" =  model_summaries,
         "contrasts_DA" = contrasts_DA,
         "contrasts_DV" = contrasts_DV,
         "discriminant_taxa_DA" = disc_vec_da,
         "discriminant_taxa_DV" = disc_vec_dv,
         "data" = data),
    class = "contrastsTest"
  )
}
