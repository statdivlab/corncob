#' Identify differentially-abundant and differentially-variable taxa
#'
#' @param formula Formula for mean, without response
#' @param phi.formula Formula for overdispersion, without response
#' @param formula_null Formula for mean under null, without response
#' @param phi.formula_null Formula for overdispersion under null, without response
#' @param data Data frame or \code{phyloseq} object
#' @param link Link function for mean, defaults to "logit"
#' @param phi.link Link function for overdispersion, defaults to "logit"
#' @param cutoff Desired type 1 error rate
#' @param fdr False discovery rate control method, defaults to "fdr"
#' @param ... Additional arguments for \code{\link{bbdml}}
#'
#' @return List of differentially-abundant taxa and differentially-variable taxa
#' @export
differentialTest <- function(formula, phi.formula,
                             formula_null, phi.formula_null,
                             data,
                             link = "logit",
                             phi.link = "logit",
                             cutoff = 0.05,
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
    colnames(out) <- c("DA","DV","Error")

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

      ### TODO
      if (ncol(X.b_i) - ncol(X.b_null_i) > 1 || ncol(X.bstar_i) - ncol(X.bstar_null_i) > 1) {
        stop("This function is in beta. \n It is currently only implemented for testing a single covariate in each link function.")
      }

      # Fit unrestricted model
      fit_unr <- try(bbdml(formula = formula_i, phi.formula = phi.formula,
                       data = data_i, link = link, phi.link = phi.link, ...),
                     silent = TRUE)
      if (class(fit_unr) == "try-error") {
        out[i, 3] <- 1
      } else {
        coef.table <- waldtest(fit_unr)
        out[i, ] <- c(coef.table[c(2, 4), 4], 0)
      } # Closes if/else for try-error check
    } ### This closes loop through taxa

    # Now have matrix of taxa with p-values and indicator for model fit
    post_fdr <- matrix(stats::p.adjust(out[,c(1,2)], method = fdr), ncol = 2)
    # Record significant taxa
    DA_vec <- rownames(out)[which(post_fdr[,1] < cutoff)]
    DV_vec <- rownames(out)[which(post_fdr[,2] < cutoff)]
    error_vec <- rownames(out)[which(post_fdr[,3] == 1)]

    return(list("p" = out, "p_fdr" = post_fdr,
                "DA" = DA_vec, "DV" = DV_vec, "error" = error_vec))

  } else { ### This closes phyloseq if
    stop("This function is in beta. \n It is currently only implemented for phyloseq objects.")
  }

}
