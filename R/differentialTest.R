#' Identify differentially-abundant and differentially-variable taxa
#'
#' @param formula Formula for mean, without response
#' @param phi.formula Formula for overdispersion, without response
#' @param formula_null Formula for mean under null, without response
#' @param phi.formula_null Formula for overdispersion under null, without response
#' @param data Data frame or matrix of count data.  Alternatively, a \code{phyloseq} object
#' @param test Character. Type of test. One of \code{"Wald"}, \code{"LRT"}
#' @param boot Boolean. Defaults to \code{FALSE}. Indicator of whether or not to use parametric bootstrap algorithm.
#' @param B Integer. Optional. Defaults to \code{NULL}. Number of bootstrap iterations. Ignored if \code{boot} is \code{FALSE}. Otherwise, if \code{NULL}, uses 1000.
#' @param sample_data Data frame or matrix. Defaults to \code{NULL}. If \code{data} is a data frame or matrix, this must be included as covariates.
#' @param taxa_are_rows Boolean. Optional. If \code{data} is a data frame or matrix, this indicates whether taxa are rows. Defaults to \code{TRUE}.
#' @param link Link function for mean, defaults to "logit"
#' @param phi.link Link function for overdispersion, defaults to "logit"
#' @param fdr_cutoff Desired type 1 error rate
#' @param fdr False discovery rate control method, defaults to "fdr"
#' @param inits Initializations for the unrestricted model to be passed to \code{\link{bbdml}}. Defaults to \code{NULL}.
#' @param inits_null Initializations for model restricted under the null to be passed to \code{\link{bbdml}}. Defaults to \code{NULL}.
#' @param ... Additional arguments for \code{\link{bbdml}}
#'
#' @details Make sure the number of columns in all of the initializations are correct! \code{inits} probably shouldn't match \code{inits_null_mu} or \code{inits_null_phi}.
#'
#' @return List of p-values and taxa names for test specified by contrast between formulas and formulas under the null.
#' @export
differentialTest <- function(formula, phi.formula,
                             formula_null, phi.formula_null,
                             data,
                             test,
                             boot,
                             B = 1000,
                             sample_data = NULL,
                             taxa_are_rows = TRUE,
                             link = "logit",
                             phi.link = "logit",
                             fdr_cutoff = 0.05,
                             fdr = "fdr",
                             inits = NULL,
                             inits_null = NULL,
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
    pvals <- rep(NA, length(taxanames))
  } else if (is.matrix(data) || is.data.frame(data)) {

    # use phyloseq
    OTU <- phyloseq::otu_table(data, taxa_are_rows = taxa_are_rows)

    # Make sample data
    sampledata <- phyloseq::sample_data(data.frame(
      sample_data,
      row.names = sample_names(OTU)
    ))

    # Make phyloseq object
    data <- phyloseq::phyloseq(OTU, sampledata)
    # Set up response
    taxanames <- phyloseq::taxa_names(data)
    pvals <- rep(NA, length(taxanames))

  } else {
    stop("Input must be either data frame, matrix, or phyloseq object!")
  }

  # check to make sure inits is of the same length
  if (!is.null(inits)) {
    ncol1 <- ncol(stats::model.matrix(object = formula, data = data.frame(sample_data(data))))
    ncol2 <- ncol(stats::model.matrix(object = phi.formula, data = data.frame(sample_data(data))))
    if (length(inits) != ncol1 + ncol2) {
      stop("inits must match number of regression parameters in formula and phi.formula!")
    }
  }
  # inits_null_mu
  if (!is.null(inits_null)) {
    ncol1 <- ncol(stats::model.matrix(object = formula_null, data = data.frame(sample_data(data))))
    ncol2 <- ncol(stats::model.matrix(object = phi.formula_null, data = data.frame(sample_data(data))))
    if (length(inits_null) != ncol1 + ncol2) {
      stop("init_null must match number of regression parameters in formula_null and phi.formula_null!")
    }
  }


    # Loop through OTU/taxa
    for (i in 1:length(taxanames)) {

      # Subset data to only select that taxa
      data_i <- convert_phylo(data, select = taxanames[i])

      # Update formula to match
      formula_i <- stats::update(formula, cbind(W, M) ~ .)
      formula_null_i <- stats::update(formula_null, cbind(W, M) ~ .)

      # Fit unrestricted model
      mod <- try(bbdml(formula = formula_i, phi.formula = phi.formula,
                    data = data_i, link = link, phi.link = phi.link,
                    inits = inits, ...), silent = TRUE)

      # Fit restricted model
      mod_null <- try(bbdml(formula = formula_null_i, phi.formula = phi.formula_null,
                       data = data_i, link = link, phi.link = phi.link,
                       inits = inits_null, ...), silent = TRUE)

      if (!("try-error" %in% c(class(mod), class(mod_null)))) {
        # If both models fit, otherwise keep as NA
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
      }
    }

    post_fdr <- stats::p.adjust(pvals, method = fdr)
    names(pvals) <- names(post_fdr) <- taxanames
    # Record significant taxa
    signif_vec <- taxanames[which(post_fdr < fdr_cutoff)]

    return(list("p" = pvals, "p_fdr" = post_fdr,
                "significant_taxa" = signif_vec))
}
