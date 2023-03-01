

#' Computes the estimated skewness of an ordinal time series
#'
#' \code{ordinal_skewness} computes the estimated skewness
#' of an ordinal time series
#'
#' @param series An OTS.
#' @param states A numerical vector containing the corresponding
#' states.
#' @param distance A function defining the underlying distance between
#' states. The Hamming, block and Euclidean distances are already
#' implemented by means of the arguments "Hamming", "Block" (default)
#' and "Euclidean". Otherwise, a function taking as input two states must
#' be provided.
#' @param normalize Logical. If \code{normalize = FALSE} (default), the value of the estimated skewness is returned. Otherwise, the function
#' returns the normalized estimated skewness.
#' @return The estimated skewness.
#' @examples
#' estimated_skewness <- ordinal_skewness(series = AustrianWages$data[[100]],
#' states = 0 : 5) # Computing the skewness estimate
#' # for one series in dataset AustrianWages using the block distance
#' @details
#' Given an OTS of length \eqn{T} with range \eqn{\mathcal{S}=\{s_0, s_1, s_2, \ldots, s_n\}} (\eqn{s_0 < s_1 < s_2 < \ldots < s_n}),
#' \eqn{\overline{X}_t=\{\overline{X}_1,\ldots, \overline{X}_T\}}, the function computes the
#' estimated skewness given by \eqn{\widehat{skew}_{d}=\sum_{i=0}^n\big(d(s_i,s_n)-d(s_i,s_0)\big)\widehat{p}_i},
#' where \eqn{d(\cdot, \cdot)} is a distance between ordinal states and \eqn{\widehat{p}_k} is the standard estimate
#' of the marginal probability for state \eqn{s_k} computed from the realization \eqn{\overline{X}_t}.
#' @encoding UTF-8
#' @author
#' Ángel López-Oriona, José A. Vilar
#' @references{
#'
#'   \insertRef{weiss2019distance}{otsfeatures}
#'
#' }
#' @export

ordinal_skewness <- function(series, states, distance = 'Block',
                              normalize = FALSE) {

  check_ots(series)
  series_length <- length(series)
  n_states <- length(states)
  distance_function <- distance

  if (distance == 'Hamming') {

    distance_function <- auxiliary_dis_hamming

  }

  if (distance == 'Block') {

    distance_function <- auxiliary_dis_block

  }

  if (distance == 'Euclidean') {

    distance_function <- auxiliary_dis_euclidean

  }

  vector_mp <- marginal_probabilities(series, states)
  d_0_m <- distance_function(states[1], states[n_states])
  vector_differences <- numeric()

  for (i in 1 : n_states) {

    vector_differences[i] <- distance_function(states[i], states[n_states]) - distance_function(states[i], states[1])

  }

  estimated_skewness <- sum(vector_mp %*% vector_differences)
  normalized_estimated_skewness <- estimated_skewness/d_0_m

  if (normalize == FALSE ){

    return(estimated_skewness)

  } else {

    return(normalized_estimated_skewness)

  }




}
