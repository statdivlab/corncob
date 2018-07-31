#' Plotting function
#'
#' @param x Output from \code{\link{bbdml}}
#' @param ... See details
#'
#' @details ... can include boolean indicator called AA for whether to plot on absolute abundance scale, which defaults to FALSE. It can also include "color", a string of a covariate in the model data, if you want the plot to be colored by color.
#'
#' @return Model plot
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#' @export
plot.bbdml <- function(x, ...) {
  input <- match.call(expand.dots = TRUE)
  mod <- x
  if (is.null(input$AA)) {
    input$AA <- FALSE
  }

  AA <- input$AA
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
  if (!AA) {
    ymin <- ymin/M
    ymax <- ymax/M
  }
  # When doing RA plot
  # If we want to create a CI for mu, uncomment (and move)

  # ginv <- switch(x$link, "logit" = invlogit)
  # g <- switch(x$link, "logit" = logit)
  #
  # # X
  # X_mu <- x$X.mu
  # # variance of beta - in Wald Test as well
  # covMat <- try(chol2inv(chol(hessian(mod))), silent = TRUE)
  # if (class(covMat) == "try-error") {
  #   stop("Singular Hessian!")
  # }
  # mu_se <- sqrt(diag(X_mu %*% covMat[1:x$np.mu,1:x$np.mu] %*% t(X_mu)))
  # ymin <- ginv(g(x$mu.resp) + qnorm(.025)*mu_se)
  # ymax <- ginv(g(x$mu.resp) + qnorm(.975)*mu_se)

  samp_names <- rownames(x$dat)


  # Fix for global bindings warnings
  index <- NULL

  if (!is.null(input$color)) {
    color <- factor(mod$dat[[input$color]])
    color_name <- input$color
    if (!AA) {
      df_RA <- data.frame(RA = W/M,
                          E_RA = mu_est,
                          SE_RA = sqrt(mu_est*(1 - mu_est)*(1 + ((M - 1)*phi_est))/M),
                          color = color,
                          #index = order(order(color)),
                          index = samp_names,
                          ymin = ymin,
                          ymax = ymax
      )
      ggplot2::ggplot(df_RA, ggplot2::aes(x = index, y = RA, color = color, group = color)) +
        ggplot2::geom_point(ggplot2::aes(color = color, group = color)) +
        #ggplot2::geom_point(ggplot2::aes(x = index, y = E_RA, color = color), pch = 2) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = ymin, ymax = ymax, color = color), width = .2,
                               position = ggplot2::position_dodge(.05)) +
        ggplot2::labs(x = "Sample", y = "Relative Abundance", title = "Model Fit", colour = color_name) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), plot.title = ggplot2::element_text(hjust = 0.5))
    } else {
      df_RA <- data.frame(RA = W,
                          E_RA = mu_est*M,
                          SE_RA = sqrt(M*mu_est*(1 - mu_est)*(1 + ((M - 1)*phi_est))),
                          color = color,
                          #index = order(order(color)),
                          index = samp_names,
                          ymin = ymin,
                          ymax = ymax
      )
      ggplot2::ggplot(df_RA, ggplot2::aes(x = index, y = RA, color = color, group = color)) +
        ggplot2::geom_point(ggplot2::aes(color = color, group = color)) +
        #ggplot2::geom_point(ggplot2::aes(x = index, y = E_RA, color = color), pch = 2) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = ymin, ymax = ymax, color = color), width = .2,
                               position = ggplot2::position_dodge(.05)) +
        ggplot2::labs(x = "Sample", y = "Absolute Abundance", title = "Model Fit", colour = color_name) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), plot.title = ggplot2::element_text(hjust = 0.5))
    }
  } else {
    if (!AA) {
      df_RA <- data.frame(RA = W/M,
                          E_RA = mu_est,
                          SE_RA = sqrt(mu_est*(1 - mu_est)*(1 + ((M - 1)*phi_est))/M),
                          #index = 1:length(M),
                          index = samp_names,
                          ymin = ymin,
                          ymax = ymax
      )

      ggplot2::ggplot(df_RA, ggplot2::aes(x = index, y = RA)) +
        ggplot2::geom_point() +
        #ggplot2::geom_point(ggplot2::aes(x = index, y = E_RA), pch = 2) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = ymin, ymax = ymax), width = .2,
                               position = ggplot2::position_dodge(.05)) +
        ggplot2::labs(x = "Sample", y = "Relative Abundance", title = "Model Fit") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), plot.title = ggplot2::element_text(hjust = 0.5))
    } else {
      df_RA <- data.frame(RA = W,
                          E_RA = mu_est*M,
                          SE_RA = sqrt(M*mu_est*(1 - mu_est)*(1 + ((M - 1)*phi_est))),
                          #index = 1:length(M),
                          index = samp_names,
                          ymin = ymin,
                          ymax = ymax
      )

      ggplot2::ggplot(df_RA, ggplot2::aes(x = index, y = RA)) +
        ggplot2::geom_point() +
        #ggplot2::geom_point(ggplot2::aes(x = index, y = E_RA), pch = 2) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = ymin, ymax = ymax), width = .2,
                               position = ggplot2::position_dodge(.05)) +
        ggplot2::labs(x = "Sample", y = "Absolute Abundance", title = "Model Fit") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), plot.title = ggplot2::element_text(hjust = 0.5))
    }
  }
}
