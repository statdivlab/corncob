#' Plotting function
#'
#' @param x Object of class \code{bbdml}.
#' @param total (Optional). Default \code{FALSE}. Boolean indicator for whether to plot on total counts scale
#' @param color (Optional). Default \code{NULL}. The sample variable to map to different colors. Can be a single character string of the variable name in \code{sample_data} or a custom supplied vector with length equal to the number of samples. Use a character vector to have \code{ggplot2} default.
#' @param shape (Optional). Default \code{NULL}. The sample variable to map to different shapes. Can be a single character string of the variable name in \code{sample_data} or a custom supplied vector with length equal to the number of samples.
#' @param facet (Optional). Default \code{NULL}. The sample variable to map to different panels in a facet grid. Must be a single character string of a variable name in \code{sample_data}.
#' @param title (Optional). Default \code{NULL}. Character string. The main title for the graphic.
#' @param B (Optional). Default \code{1000}. Integer. Number of bootstrap simulations for prediction intervals. Use \code{B = 0} for no prediction intervals.
#' @param sample_names (Optional). Default \code{TRUE}. Boolean. If \code{FALSE}, remove sample names from the plot.
#' @param data_only (Optional). Default \code{FALSE}. Boolean. If \code{TRUE}, only returns data frame.
#' @param ... There are no optional parameters at this time.
#'
#' @return Object of class \code{ggplot}. Plot of \code{bbdml} model fit with 95% prediction intervals.
#'
#' @examples
#' \donttest{
# data(soil_phylo)
# soil <- soil_phylo %>%
# phyloseq::subset_samples(DayAmdmt %in% c(11,21)) %>%
# phyloseq::tax_glom("Phylum")
# mod <- bbdml(formula = OTU.1 ~ DayAmdmt,
# phi.formula = ~ DayAmdmt,
# data = soil)
# plot(mod, color = "DayAmdmt")
#' }
#' @export
plot.bbdml <- function(x, total = FALSE, color = NULL, shape = NULL, facet = NULL, title = NULL, B = 1000, sample_names = TRUE, data_only = FALSE, ...) {
  # input <- match.call(expand.dots = TRUE)
  mod <- x

  M <- mod$M
  W <- mod$W


  if (B != 0) {
    ymin <- ymax <- rep(NA, length(M))


    sims <- matrix(NA, nrow = B, ncol = length(W))
    newdat <- mod$dat
    for (i in 1:B) {
      sim <- simulate(mod, nsim = length(W))
      newdat$W <- sim
      refit <- suppressWarnings(bbdml(mod$formula, phi.formula = mod$phi.formula,
                                      link = mod$link, phi.link = mod$phi.link,
                                      inits = mod$inits,
                                      data = newdat))
      sims[i,] <- simulate(refit, nsim = length(W))
    }

    predint <- apply(sims, 2, stats::quantile, c(.025, .975))
    ymin <- predint[1,]
    ymax <- predint[2,]
  } else {
    ymin <- ymax <- W
  }

  resp <- W
  if (!total) {
    ymin <- ymin/M
    ymax <- ymax/M
    resp <- W/M
  }

  samp_names <- names(W)
  dat_noNA <- mod$dat[samp_names,]

  df <- data.frame(RA = resp,
                      samples = samp_names,
                      ymin = ymin,
                      ymax = ymax
  )

  if (!(data_only)) {

  my_ord_str <- ""
  custom_color <- custom_shape <- FALSE
  color_name <- shape_name <- NULL
  if (!is.null(color)) {
    if (length(color) == 1) {
      df[[color]] <- factor(dat_noNA[[color]])
      color_name <- color
      my_ord_str <- paste(my_ord_str, df[[color]], sep = "_")
    } else if (length(color) == nrow(df)) {
      df[["color"]] <- color
      color <- color_name <- "color"
      custom_color <- TRUE
    } else {
      stop("color must either match a variable or be a custom vector of correct length!")
    }
    colvar <- utils::tail(colnames(df), 1)
  } else {
    df[["color"]] <- NA
    colvar <- "color"
  } # End if (!is.null(color))

  if (!is.null(shape)) {
    if (length(shape) == 1) {
      df[[shape]] <- factor(dat_noNA[[shape]])
      shape_name <- shape
      my_ord_str <- paste(my_ord_str, df[[shape]], sep = "_")
    } else if (length(shape) == nrow(df)) {
      df[["shape"]] <- shape
      shape <- shape_name <- "shape"
      custom_shape <- TRUE
    } else {
      stop("shape must either match a variable or be a custom vector of correct length!")
    }
    shapevar <- utils::tail(colnames(df), 1)
  } else {
    df[["shape"]] <- NA
    shapevar <- "shape"
  } # End if (!is.null(shape))

  if (!is.null(facet)) {
    df[[facet]] <- factor(dat_noNA[[facet]])
  }

  # reorder
  #my_ord_str <- paste(my_ord_str, df$samples, sep = "_")
  df$order <- factor(df$samples, levels = dplyr::arrange(df, df[[colvar]],
                                                  df[[shapevar]], samples)$samples)

  ylab_tmp <- ifelse(!total, "Relative Abundance", "Total Counts")

  aes_map <- ggplot2::aes_string(x = "order", y = "RA", colour = color, shape = shape, labs = "samples")
  my_gg <- ggplot2::ggplot(df, aes_map) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = ymin, ymax = ymax), width = .2) +
    ggplot2::labs(title = title, x = "", y = ylab_tmp, colour = color_name, shape = shape_name) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (custom_color) {
    my_gg <- my_gg + ggplot2::guides(colour = FALSE)
  }
  if (custom_shape) {
    my_gg <- my_gg + ggplot2::scale_shape_identity()
  }


  if (!is.null(facet)) {
    my_gg <- my_gg + ggplot2::facet_grid(paste0("~", facet), scales = "free_x", space = "free_x", labeller = ggplot2::label_both)
  }

  if (!(sample_names)) {
    my_gg <- my_gg + ggplot2::theme(axis.text.x = ggplot2::element_blank())
  }
  my_gg
  } else {
    return(df)
  }
}

# remove global variable NOTE
utils::globalVariables("samples", add = FALSE)
