

#' Computes the conditional probabilities of an ordinal time series
#'
#' \code{conditional_probabilities} returns a matrix with the conditional
#' probabilities of an ordinal time series
#'
#' @param series An OTS.
#' @param lag The considered lag (default is 1).
#' @param states A numerical vector containing the corresponding
#' states.
#' @return A matrix with the conditional probabilities.
#' @examples
#' matrix_cp <- conditional_probabilities(series = AustrianWages$data[[100]],
#' states = 0 : 5) # Computing the matrix of
#' # conditional probabilities for one series in dataset AustrianWages
#' @details
#' Given an OTS of length \eqn{T} with range \eqn{\mathcal{S}=\{s_0, s_1, s_2, \ldots, s_n\}} (\eqn{s_0 < s_1 < s_2 < \ldots < s_n}),
#' \eqn{\overline{X}_t=\{\overline{X}_1,\ldots, \overline{X}_T\}}, the function computes the
#' matrix \eqn{\widehat{\boldsymbol P}^c(l) = \big(\widehat{p}^c_{i-1j-1}(l)\big)_{1 \le i, j \le n+1}},
#' with \eqn{\widehat{p}^c_{ij}(l)=\frac{TN_{ij}(l)}{(T-l)N_i}}, where
#' \eqn{N_i} is the number of elements equal to \eqn{s_i} in the realization \eqn{\overline{X}_t} and \eqn{N_{ij}(l)} is the number
#' of pairs \eqn{(\overline{X}_t, \overline{X}_{t-l})=(s_i,s_j)} in the realization \eqn{\overline{X}_t}.
#' @encoding UTF-8
#' @author
#' Ángel López-Oriona, José A. Vilar
#' @references{
#'
#'   \insertRef{weiss2019distance}{otsfeatures}
#'
#' }
#' @export

conditional_probabilities <- function(series, lag = 1, states) {

  check_ots(series)
  n_states <- length(states)
  matrix_probs <- matrix(0, n_states, n_states)

  vector_mp <- marginal_probabilities(series = series, states = states)
  matrix_jp <- joint_probabilities(series = series, lag = lag,
                                   states = states)
  matrix_mp <- base::matrix(vector_mp, nrow = n_states, ncol = n_states, byrow = TRUE)
  matrix_conditional_probabilities <- matrix_jp/matrix_mp




  return(matrix_conditional_probabilities)

}
