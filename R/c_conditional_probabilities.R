

#' Computes the cumulative conditional probabilities of an ordinal time series
#'
#' \code{c_conditional_probabilities} returns a matrix with the cumulative conditional
#' probabilities of an ordinal time series
#'
#' @param series An OTS.
#' @param lag The considered lag (default is 1).
#' @param states A numerical vector containing the corresponding
#' states.
#' @return A matrix with the conditional probabilities.
#' @examples
#' matrix_ccp <- c_conditional_probabilities(series = AustrianWages$data[[100]],
#' states = 0 : 5) # Computing the matrix of
#' # cumulative conditional probabilities for one series in dataset AustrianWages
#' @details
#' Given an OTS of length \eqn{T} with range \eqn{\mathcal{S}=\{s_0, s_1, s_2, \ldots, s_n\}} (\eqn{s_0 < s_1 < s_2 < \ldots < s_n}),
#' \eqn{\overline{X}_t=\{\overline{X}_1,\ldots, \overline{X}_T\}}, the function computes the
#' matrix \eqn{\widehat{\boldsymbol F}^c(l) = \big(\widehat{f}^c_{i-1j-1}(l)\big)_{1 \le i, j \le n}},
#' with \eqn{\widehat{f}^c_{ij}(l)=\frac{TN_{ij}(l)}{(T-l)N_i}}, where
#' \eqn{N_i} is the number of elements less one or equal to \eqn{s_i} in the realization \eqn{\overline{X}_t} and \eqn{N_{ij}(l)} is the number
#' of pairs \eqn{(\overline{X}_t, \overline{X}_{t-l})} in the realization \eqn{\overline{X}_t}
#' such that \eqn{\overline{X}_t \le s_i} and \eqn{\overline{X}_{t-l} \le s_j}.
#' @encoding UTF-8
#' @author
#' Ángel López-Oriona, José A. Vilar
#' @references{
#'
#'   \insertRef{weiss2019distance}{otsfeatures}
#'
#' }
#' @export

c_conditional_probabilities <- function(series, lag = 1, states) {

  check_ots(series)
  n_states <- length(states)
  matrix_probs <- matrix(0, n_states, n_states)

  vector_mp <- c_marginal_probabilities(series = series, states = states)
  matrix_jp <- c_joint_probabilities(series = series, lag = lag,
                                   states = states)
  matrix_mp <- base::matrix(vector_mp, nrow = n_states - 1, ncol = n_states - 1, byrow = TRUE)
  matrix_conditional_probabilities <- matrix_jp/matrix_mp




  return(matrix_conditional_probabilities[-n_states, -n_states])

}
