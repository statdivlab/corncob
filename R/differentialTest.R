#' Identify differentially-abundant and differentially-variable taxa
#'
#' @param formula Formula for mean, without response
#' @param phi.formula Formula for overdispersion, without response
#' @param formula_null Formula for mean under null, without response
#' @param phi.formula_null Formula for overdispersion under null, without response
#' @param data Data frame, matrix, or \code{phyloseq} object
#' @param link Link function for mean, defaults to "logit"
#' @param phi.link Link function for overdispersion, defaults to "logit"
#' @param fdr_cutoff Desired type 1 error rate
#' @param fdr False discovery rate control method, defaults to "fdr"
#' @param ... Additional arguments for \code{\link{bbdml}}
#'
#' @details Note that if you are testing a single covariate, this function will use a Wald test. If you are testing multiple covariates, this function will use a likelihood ratio test.
#'
#' @return List of differentially-abundant taxa and differentially-variable taxa
#' @export
differentialTest <- function(formula, phi.formula,
                             formula_null, phi.formula_null,
                             data,
                             link = "logit",
                             phi.link = "logit",
                             fdr_cutoff = 0.05,
                             fdr = "fdr",
                             ...) {
  # Record call
  call <- match.call(expand.dots = FALSE)
  # Record mu link
  link <- match.arg(link, choices = "logit")
  # Record phi link
  phi.link <- match.arg(phi.link, choices = c("fishZ", "logit"))

  # Convert phyloseq objects
  if ("phyloseq" %in% class(data)) {

    # Set up response
    out <- matrix(NA, ncol = 3, nrow = nrow(phyloseq::otu_table(data)))
    rownames(out) <- rownames(phyloseq::otu_table(data))
  } else if (is.matrix(data) || is.data.frame(data)) {
    out <- matrix(NA, ncol = 3, nrow = nrow(data))
    rownames(out) <- rownames(data)
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
                             data = data_i, link = link, phi.link = phi.link, ...),
                       silent = TRUE)

        # If doesn't load, don't need to test other model
        if (class(fit_unr) == "try-error") {
          out[i, 3] <- 1
        } else {
          # Fit restricted model  for mu
          fit_res_mu <- try(bbdml(formula = formula_null_i, phi.formula = phi.formula,
                               data = data_i, link = link, phi.link = phi.link, ...),
                         silent = TRUE)
          p.val_mu <- try(lrtest(fit_res_mu, fit_unr), silent = TRUE)

          # Fit restricted model for phi
          fit_res_phi <- try(bbdml(formula = formula_i, phi.formula = phi.formula_null,
                                  data = data_i, link = link, phi.link = phi.link, ...),
                            silent = TRUE)
          p.val_phi <- try(lrtest(fit_res_phi, fit_unr), silent = TRUE)

          # Begin testing fit_res_mu
          if (class(fit_res_mu) == "try-error" || class(p.val_mu) == "try-error") {
            out[i, 3] <- 1
          } else {
            out[i, 1] <- lrtest(fit_res_mu, fit_unr)
          } # if fit_res_mu or p.val breaks

          # Begin testing fit_res_phi
          if (class(fit_res_phi) == "try-error" || class(p.val_phi) == "try-error") {
            out[i, 3] <- 1
          } else {
            out[i, 2] <- lrtest(fit_res_phi, fit_unr)
          } # if fit_res_phi or p.val breaks

          # Put a 0 if neither breaks
          if (is.na(out[i, 3])) {
            out[i, 3] <- 0
          } # end if need to replace warning with 0

        } # end else after testing that fit_unr fit

      } # End multiple testing if

      # Fit unrestricted model
      fit_unr <- try(bbdml(formula = formula_i, phi.formula = phi.formula,
                       data = data_i, link = link, phi.link = phi.link, ...),
                     silent = TRUE)

      coef.table <- try(waldtest(fit_unr), silent = TRUE)

      if (class(fit_unr) == "try-error" || class(coef.table) == "try-error") {
        out[i, 3] <- 1
      } else {
          out[i, ] <- c(coef.table[c(2, 4), 4], 0)
      } # Closes if/else for try-error check
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
