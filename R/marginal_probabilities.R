

#' Computes the marginal probabilities of an ordinal time series
#'
#' \code{marginal_probabilities} returns a vector with the marginal
#' probabilities of an ordinal time series
#'
#' @param series An OTS (numerical vector with integers).
#' @param states A numerical vector containing the corresponding
#' states
#' @return A vector with the marginal probabilities.
#' @examples
#' vector_mp <- marginal_probabilities(series = AustrianWages$data[[100]],
#' states = 0 : 5) # Computing the vector of
#' # marginal probabilities for one series in dataset AustrianWages
#' @details
#' Given an OTS of length \eqn{T} with range \eqn{\mathcal{S}=\{s_0, s_1, s_2, \ldots, s_n\}} (\eqn{s_0 < s_1 < s_2 < \ldots < s_n}),
#' \eqn{\overline{X}_t=\{\overline{X}_1,\ldots, \overline{X}_T\}}, the function computes the
#' vector \eqn{\widehat{\boldsymbol p} =(\widehat{p}_0, \ldots, \widehat{p}_n)},
#' with \eqn{\widehat{p}_i=\frac{N_i}{T}}, where \eqn{N_i} is the number
#' of elements equal to \eqn{s_i} in the realization \eqn{\overline{X}_t}.
#' @encoding UTF-8
#' @author
#' Ángel López-Oriona, José A. Vilar
#' @references{
#'
#'   \insertRef{weiss2019distance}{otsfeatures}
#'
#' }
#' @export

marginal_probabilities <- function(series, states) {

  check_ots(series)
  series_length <- length(series) # Series length
  n_states <- length(states) # Number of states in the dataset


  # Computing the marginal probabilities in the series

  marginal_probabilities <- numeric()

  for (i in 1 : n_states) {

    count_i <- sum(series == states[[i]])
    marginal_probabilities[i] <- count_i/series_length

  }

  return(marginal_probabilities)

}
