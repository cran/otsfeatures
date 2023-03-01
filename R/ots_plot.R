

#' Constructs an ordinal time series plot
#'
#' \code{ots_plot} constructs an ordinal time series plot
#'
#' @param series An OTS.
#' @param states A numerical vector containing the corresponding
#' states.
#' @param title The title of the graph.
#' @param labels The labels of the graph.
#' @return The ordinal time series plot.
#' @examples
#' ordinal_time_series_plot <- ots_plot(series = AustrianWages$data[[100]],,
#' states = 0 : 5) # Constructs an ordinal
#' # time series plot for one series in
#' # dataset AustrianWages
#' @details
#' Constructs an ordinal time series plot for a given OTS.
#' @encoding UTF-8
#' @author
#' Ángel López-Oriona, José A. Vilar
#' @references{
#'
#'   \insertRef{weiss2018introduction}{otsfeatures}
#'
#' }
#' @export

ots_plot <- function(series, states, title = 'Time series plot', labels = NULL) {

  x <- y <- NULL
  check_ots(series)
  series_length <- length(series)
  n_states <- length(states)
  numeric_series <- numeric(series_length)

  for (i in 1 : n_states) {

    indexes_i <- which(series == states[i])
    numeric_series[indexes_i] <- i

  }

  if (is.null(labels)) {

  labels <- numeric(n_states)



  for (i in 1 : n_states) {

  labels[i] <- paste0('s', i-1)

  }

  }

  df_plot <- data.frame(x = 1 : series_length, y = numeric_series)

  plot_ots <- ggplot2::ggplot(df_plot, ggplot2::aes(x = x, y = y)) + ggplot2::geom_line(size = 1, col = 'blue') +
    ggplot2::geom_point(size = 1.5, col = 'blue') + ggplot2::xlab('Time') +
    ggplot2::scale_y_continuous(breaks = as.numeric(states + 1), labels = labels) +
    ggplot2::ylab('') + ggplot2::theme(axis.text = ggplot2::element_text(size = 11),
                                       axis.title.x = ggplot2::element_text(size = 12),
                                       plot.title = ggplot2::element_text(hjust = 0.5, size = 12)) +
    ggplot2::ggtitle(title)

  return(plot_ots)

}

