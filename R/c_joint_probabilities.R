

#' Computes the cumulative joint probabilities of an ordinal time series
#'
#' \code{c_joint_probabilities} returns a matrix with the cumulative joint
#' probabilities of an ordinal time series
#'
#' @param series An OTS.
#' @param lag The considered lag (default is 1).
#' @param states A numerical vector containing the corresponding
#' states.
#' @return A matrix with the jcumulative oint probabilities.
#' @examples
#' matrix_cjp <- c_joint_probabilities(series = AustrianWages$data[[100]],
#' states = 0 : 5) # Computing the matrix of
#' # cumulative joint probabilities for one series in dataset AustrianWages
#' @details
#' Given an OTS of length \eqn{T} with range \eqn{\mathcal{S}=\{s_0, s_1, s_2, \ldots, s_n\}} (\eqn{s_0 < s_1 < s_2 < \ldots < s_n}),
#' \eqn{\overline{X}_t=\{\overline{X}_1,\ldots, \overline{X}_T\}}, the function computes the
#' matrix \eqn{\widehat{\boldsymbol F}(l) = \big(\widehat{f}_{i-1j-1}(l)\big)_{1 \le i, j \le n}},
#' with \eqn{\widehat{f}_{ij}(l)=\frac{N_{ij}(l)}{T-l}}, where \eqn{N_{ij}(l)} is the number
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

c_joint_probabilities <- function(series, lag = 1, states) {

  check_ots(series)
  n_states <- length(states)
  matrix_joint_probabilities <- base::matrix(0, n_states, n_states)

  for (i in 1 : n_states) {

    for (j in 1 : n_states) {

      matrix_joint_probabilities[i, j] <- f_i_j_k_function(series = series, i_stat  = states[i],
                                                           j_stat = states[j], k = lag)

    }

  }

  return(matrix_joint_probabilities[-n_states, -n_states])

}
