#' Plotting function (in beta)
#'
#' @param x Output from \code{\link{bbdml}}
#' @param ... Optional boolean indicator called RA for whether to plot on relative abundance scale. Defaults to FALSE.
#'
#' @return Model plot
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#' @export
plot.bbdml <- function(x, ...) {
  mod <- x
  if (!exists("RA")) {
    RA <- FALSE
  }
  mu_est <- mod$mu.resp
  phi_est <- mod$phi.resp
  M <- mod$M
  W <- mod$W
  if (RA) {
    df_RA <- data.frame(RA = W/M,
                        E_RA = mu_est,
                        SE_RA = sqrt(mu_est*(1 - mu_est)*(1 + ((M - 1)*phi_est))/M),
                        index = 1:length(M)
    )

    ggplot2::ggplot(df_RA, ggplot2::aes(x = index, y = RA)) +
      ggplot2::geom_point() +
      ggplot2::geom_point(ggplot2::aes(x = index, y = E_RA), pch = 2) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = E_RA - 2*SE_RA, ymax = E_RA + 2*SE_RA), width = .2,
                    position = ggplot2::position_dodge(.05)) +
      ggplot2::labs(x = "Sample", y = "Relative Abundance", title = "BBD Model Fit") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  } else {
    df_RA <- data.frame(RA = W,
                        E_RA = mu_est*M,
                        SE_RA = sqrt(M*mu_est*(1 - mu_est)*(1 + ((M - 1)*phi_est))),
                        index = 1:length(M)
    )

    ggplot2::ggplot(df_RA, ggplot2::aes(x = index, y = RA)) +
      ggplot2::geom_point() +
      ggplot2::geom_point(ggplot2::aes(x = index, y = E_RA), pch = 2) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = E_RA - 2*SE_RA, ymax = E_RA + 2*SE_RA), width = .2,
                    position = ggplot2::position_dodge(.05)) +
      ggplot2::labs(x = "Sample", y = "Absolute Abundance", title = "BBD Model Fit") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }


}
