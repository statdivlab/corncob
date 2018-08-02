#' Identify differentially-abundant and differentially-variable taxa
#'
#' @param formula Formula for mean, without response
#' @param phi.formula Formula for overdispersion, without response
#' @param formula_null Formula for mean under null, without response
#' @param phi.formula_null Formula for overdispersion under null, without response
#' @param data Data frame or matrix of count data.  Alternatively, a \code{phyloseq} object
#' @param sample_data Data frame or matrix. Defaults to \code{NULL}. If \code{data} is a data frame or matrix, this must be included as covariates.
#' @param taxa_are_rows Boolean. Optional. If \code{data} is a data frame or matrix, this indicates whether taxa are rows. Defaults to \code{TRUE}.
#' @param link Link function for mean, defaults to "logit"
#' @param phi.link Link function for overdispersion, defaults to "logit"
#' @param fdr_cutoff Desired type 1 error rate
#' @param fdr False discovery rate control method, defaults to "fdr"
#' @param inits Initializations for the unrestricted model to be passed to \code{\link{bbdml}}. Defaults to \code{NULL}.
#' @param inits_null_mu Initializations for model restricted under \code{formula_null} to be passed to \code{\link{bbdml}}. Defaults to \code{NULL}.
#' @param inits_null_phi Initializations for model restricted under \code{phi.formula_null} to be passed to \code{\link{bbdml}}. Defaults to \code{NULL}.
#' @param ... Additional arguments for \code{\link{bbdml}}
#'
#' @details Note that if you are testing a single covariate, this function will use a Wald test. If you are testing multiple covariates, this function will use a likelihood ratio test. If you are testing multiple variables. Make sure the number of columns in all of the initializations are correct! \code{inits} probably shouldn't match \code{inits_null_mu} or \code{inits_null_phi}.
#'
#' @return List of differentially-abundant taxa and differentially-variable taxa
#' @export
differentialTest <- function(formula, phi.formula,
                             formula_null, phi.formula_null,
                             data, sample_data = NULL,
                             taxa_are_rows = TRUE,
                             link = "logit",
                             phi.link = "logit",
                             fdr_cutoff = 0.05,
                             fdr = "fdr",
                             inits = NULL,
                             inits_null_mu = NULL,
                             inits_null_phi = NULL,
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
    out <- matrix(NA, ncol = 3, nrow = length(phyloseq::taxa_names(data)))
    rownames(out) <- rownames(phyloseq::otu_table(data))
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
    out <- matrix(NA, ncol = 3, nrow = length(phyloseq::taxa_names(data)))
    rownames(out) <- phyloseq::taxa_names(data)
  } else {
    stop("Input must be either data frame, matrix, or phyloseq object!")
  }

  colnames(out) <- c("DA","DV","Warning")

    # Loop through OTU/taxa
    for (i in 1:nrow(out)) {

      # Subset data to only select that taxa
      data_i <- convert_phylo(data, select = rownames(out)[i])

      # Update formula to match
      formula_i <- stats::update(formula, cbind(W, M) ~ .)
      formula_null_i <- stats::update(formula_null, cbind(W, M) ~ .)

      # mu model matrix
      X.b_i <- stats::model.matrix(object = formula_i, data = data_i)
      # phi model matrix
      X.bstar_i <- stats::model.matrix(object = phi.formula, data = data_i)
      # mu model matrix under null
      X.b_null_i <- stats::model.matrix(object = formula_null_i, data = data_i)
      # phi model matrix
      X.bstar_null_i <- stats::model.matrix(object = phi.formula_null, data = data_i)

      if (ncol(X.b_i) - ncol(X.b_null_i) > 1 || ncol(X.bstar_i) - ncol(X.bstar_null_i) > 1) {
        # Begin multiple testing if

        # Fit unrestricted model
        fit_unr <- try(bbdml(formula = formula_i, phi.formula = phi.formula,
                             data = data_i, link = link, phi.link = phi.link,
                             inits = inits, ...),
                       silent = TRUE)

        # If doesn't load, don't need to test other model
        if (class(fit_unr) == "try-error") {
          out[i, 3] <- 1
        } else {
          # Fit restricted model  for mu
          fit_res_mu <- try(bbdml(formula = formula_null_i, phi.formula = phi.formula,
                               data = data_i, link = link, phi.link = phi.link,
                               inits = inits_null_mu, ...),
                         silent = TRUE)
          p.val_mu <- try(lrtest(fit_res_mu, fit_unr), silent = TRUE)

          # Fit restricted model for phi
          fit_res_phi <- try(bbdml(formula = formula_i, phi.formula = phi.formula_null,
                                  data = data_i, link = link, phi.link = phi.link,
                                  inits = inits_null_phi, ...),
                            silent = TRUE)
          p.val_phi <- try(lrtest(fit_res_phi, fit_unr), silent = TRUE)

          # Begin testing  for try_error
          if ("try-error" %in% c(class(fit_res_mu), class(p.val_mu), class(fit_res_phi), class(p.val_phi))) {
            out[i, 3] <- 1
          } else {
            out[i, 3] <- 0
          } # if fit_res_mu or p.val breaks

          if (class(p.val_mu) != "try-error") {
            out[i, 1] <- p.val_mu
          }
          if (class(p.val_phi) != "try-error") {
            out[i, 2] <- p.val_phi
          }
        } # end else after testing that fit_unr fit
        # End multiple covariates LRT if
      } else {

        # Fit unrestricted model
        fit_unr <- try(bbdml(formula = formula_i, phi.formula = phi.formula,
                         data = data_i, link = link, phi.link = phi.link,
                         inits = inits, ...),
                       silent = TRUE)

        coef.table <- try(waldtest(fit_unr), silent = TRUE)

        if (class(fit_unr) == "try-error" || class(coef.table) == "try-error") {
          out[i, 3] <- 1
        } else {
            out[i, ] <- c(coef.table[c(2, 4), 4], 0)
        } # Closes if/else for try-error check
      } # End wald test
    } ### This closes loop through taxa

    # Now have matrix of taxa with p-values and indicator for model fit
    post_fdr <- matrix(stats::p.adjust(out[,c(1,2)], method = fdr), ncol = 2)
    colnames(post_fdr) <- colnames(out)[1:2]
    rownames(post_fdr) <- rownames(out)
    # Record significant taxa
    DA_vec <- rownames(out)[which(post_fdr[,1] < fdr_cutoff)]
    DV_vec <- rownames(out)[which(post_fdr[,2] < fdr_cutoff)]
    warning_vec <- rownames(out)[which(out[,3] == 1)]

    return(list("p" = out, "p_fdr" = post_fdr,
                "DA" = DA_vec, "DV" = DV_vec, "warning" = warning_vec))
}
