

#' Constructs a serial dependence plot based on the ordinal Cohen's kappa
#' considering the block distance
#'
#' \code{plot_ordinal_cohens_kappa} constructs a serial dependence plot of an ordinal
#' time series based on the ordinal Cohen's kappa considering the block distance
#'
#' @param series An OTS.
#' @param states A numerical vector containing the corresponding
#' states.
#' @param max_lag The maximum lag represented in the plot (default is 10).
#' @param alpha The significance level for the corresponding hypothesis test (default is 0.05).
#' @param plot Logical. If \code{plot = TRUE} (default), returns the serial dependence
#' plot. Otherwise, returns a list with the values of the ordinal Cohens's kappa, the critical
#' value and the corresponding p-values.
#' @param title The title of the graph.
#' @param bar_width The width of the corresponding bars.
#' @param ... Additional parameters for the function.
#' @return If \code{plot = TRUE} (default), returns the serial dependence plot based on the ordinal Cohens's kappa. Otherwise, the function
#' returns a list with the values of the ordinal Cohens's kappa, the critical
#' value and the corresponding p-values.
#' @examples
#' plot_ock <- plot_ordinal_cohens_kappa(series = AustrianWages$data[[100]],
#' states = 0 : 5, max_lag = 3) # Representing
#' # the serial dependence plot
#' list_ck <- plot_ordinal_cohens_kappa(series = AustrianWages$data[[100]],
#' states = 0 : 5, max_lag = 3, plot = FALSE) # Obtaining
#' # the values of the ordinal Cohens's kappa, the critical value and the p-values
#' @details
#' Constructs a serial dependence plot based on the ordinal Cohens's kappa, \eqn{\widehat{\kappa}_d(l)},
#' for several lags, where \eqn{d} is the block distance between ordinal states, that is, \eqn{d(s_i, s_j)=|i-j|} for two states \eqn{s_i} and \eqn{s_j}.
#' A dashed lined is incorporated indicating the critical value
#' of the test based on the following asymptotic approximation (under the i.i.d. assumption):
#' \deqn{\sqrt{\frac{T\widehat{disp}_d^2}{4\sum_{k,l=0}^{n-1}(\widehat{f}_{min\{k,l\}}-\widehat{f}_k\widehat{f}_l)^2}}\bigg(\widehat{\kappa}_d(l)+\frac{1}{T}\bigg)\sim N\big(0, 1\big),} where \eqn{T} is the series length,
#' \eqn{\widehat{f_k}} is the estimated cumulative probability for state \eqn{s_k}
#' and \eqn{\widehat{disp}_d} is the DIVC estimate of the dispersion.
#' @encoding UTF-8
#' @author
#' Ángel López-Oriona, José A. Vilar
#' @references{
#'
#'   \insertRef{weiss2019distance}{otsfeatures}
#'
#' }
#' @export

plot_ordinal_cohens_kappa <- function(series, states, max_lag = 10, alpha = 0.05, plot = TRUE,
                              title = 'Serial dependence plot', bar_width = 0.12,...) {

  x <- y <- NULL
  check_ots(series)
  auxiliary_list <- auxiliary_ordinal_cohens_kappa_do1(series, states,
                                           max_lag , alpha = alpha)

  df_plot_1 <- data.frame(x = 1 : max_lag, y = auxiliary_list$values_ordinal_cohens_kappa)
  df_plot_2 <- data.frame(x = 1 : (max_lag), y = rep(auxiliary_list$critical_value, max_lag))

  plot_ordinal_cohens_kappa <- ggplot2::ggplot(data = df_plot_1, mapping = ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_bar(stat = "identity", position = "identity", width = bar_width, fill = 'orange') +
    ggplot2::scale_x_continuous(breaks = 1 : max_lag) +
    ggplot2::ggtitle(title) +
    ggplot2::geom_line(data = df_plot_2, mapping = ggplot2::aes(x = x, y = y), linetype = 2, size = 0.7) +
    ggplot2::geom_line(data = df_plot_2, mapping = ggplot2::aes(x = x, y = -y), linetype = 2, size = 0.7) +
    ggplot2::xlab('Lag') + ggplot2::ylab(latex2exp::TeX("Ordinal Cohen's $\\kappa$")) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12),
                   axis.text.y = ggplot2::element_text(size = 11),
                   axis.title = ggplot2::element_text(size = 12),
                   plot.title = ggplot2::element_text(hjust = 0.5, size = 12),...)

  if (plot == TRUE) {

    return(plot_ordinal_cohens_kappa)

  } else {

    return_list <- list(values = auxiliary_list$values_ordinal_cohens_kappa,
                        p_values = auxiliary_list$vector_p_values,
                        critical_values = auxiliary_list$critical_value)

    return(return_list)

  }



}
