

#' Computes the  estimated index of ordinal variation (IOV) of an ordinal time series
#'
#' \code{index_ordinal_variation} computes the estimated index of ordinal variation
#' of an ordinal time series
#'
#' @param series An OTS.
#' @param states A numerical vector containing the corresponding
#' states.
#' @return The estimated IOV.
#' @examples
#' estimated_iov <- index_ordinal_variation(series = AustrianWages$data[[100]],
#' states = 0 : 5) # Computing the estimate of the IOV
#' # for one series in dataset AustrianWages
#' @details
#' Given an OTS of length \eqn{T} with range \eqn{\mathcal{S}=\{s_0, s_1, s_2, \ldots, s_n\}} (\eqn{s_0 < s_1 < s_2 < \ldots < s_n}),
#' \eqn{\overline{X}_t=\{\overline{X}_1,\ldots, \overline{X}_T\}}, the function computes the
#' estimated IOV given by \eqn{\widehat{IOV}=\frac{4}{n}\sum_{k=1}^{n-1}\widehat{f}_k(1-\widehat{f}_k)},
#' where \eqn{\widehat{f}_k} is the standard estimate of the cumulative marginal probability
#' for state \eqn{s_k} computed from the series \eqn{\overline{X}_t}.
#' @encoding UTF-8
#' @author
#' Ángel López-Oriona, José A. Vilar
#' @references{
#'
#'   \insertRef{weiss2019distance}{otsfeatures}
#'
#' }
#' @export

index_ordinal_variation <- function(series, states) {

  check_ots(series)
  n_states <- length(states)
  c_probabilities <- c_marginal_probabilities(series, states)
  v_c_probabilities <- c_probabilities * (1 - c_probabilities)

  return(4/(n_states - 1) * sum(v_c_probabilities))

}
