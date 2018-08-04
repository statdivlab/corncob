#' Plotting function
#'
#' @param x Output from \code{\link{bbdml}} of class \code{bbdml}.
#' @param AA (Optional). Default \code{FALSE}. Boolean indicator for whether to plot on absolute abundance scale
#' @param color (Optional). Default \code{NULL}. The sample variable to map to different colors. Can be a single character string of the variable name in \code{sample_data} or a custom supplied vector with length equal to the number of samples. Use a character vector to have \code{ggplot2} default.
#' @param shape (Optional). Default \code{NULL}. The sample variable to map to different shapes. Can be a single character string of the variable name in \code{sample_data} or a custom supplied vector with length equal to the number of samples.
#' @param facet (Optional). Default \code{NULL}. The sample variable to map to different panels in a facet grid. Must be a single character string of a variable name in \code{sample_data}.
#' @param title (Optional). Default NULL. Character string. The main title for the graphic.
#' @param ... There are no optional parameters at this time.
#'
#' @return Model plot
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#' @export
plot.bbdml <- function(x, AA = FALSE, color = NULL, shape = NULL, facet = NULL, title = NULL, ...) {
  # input <- match.call(expand.dots = TRUE)
  mod <- x

  mu_est <- mod$mu.resp
  phi_est <- mod$phi.resp
  M <- mod$M
  W <- mod$W

  ymin <- ymax <- rep(NA, length(M))

  for (i in 1:length(M)) {
    HPD <- HPDbetabinom(percent = 0.95, size = M[i], mu = mu_est[i], phi = phi_est[i])
    ymin[i] <- HPD$lower
    ymax[i] <- HPD$upper
  }
  resp <- W
  if (!AA) {
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
  } # End if (!is.null(shape))

  if (!is.null(facet)) {
    df[[facet]] <- factor(dat_noNA[[facet]])
  }

  # reorder
  my_ord_str <- paste(my_ord_str, df$samples, sep = "_")
  df$order <- factor(df$samples, levels = df$samples[order(my_ord_str)])

  ylab_tmp <- ifelse(!AA, "Relative Abundance", "Absolute Abundance")

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
  my_gg
}
