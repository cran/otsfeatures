

#' Computes the estimated asymmetry of an ordinal time series
#'
#' \code{ordinal_asymmetry} computes the estimated asymmetry
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
#' @param normalize Logical. If \code{normalize = FALSE} (default), the value of the estimated asymmetry is returned. Otherwise, the function
#' returns the normalized estimated asymmetry.
#' @return The estimated asymmetry.
#' @examples
#' estimated_asymmetry <- ordinal_asymmetry(series = AustrianWages$data[[100]],
#' states = 0 : 5) # Computing the asymmetry estimate
#' # for one series in dataset AustrianWages using the block distance
#' @details
#' Given an OTS of length \eqn{T} with range \eqn{\mathcal{S}=\{s_0, s_1, s_2, \ldots, s_n\}} (\eqn{s_0 < s_1 < s_2 < \ldots < s_n}),
#' \eqn{\overline{X}_t=\{\overline{X}_1,\ldots, \overline{X}_T\}}, the function computes the
#' estimated asymmetry given by \eqn{\widehat{asym}_{d}=\widehat{\boldsymbol p}^\top (\boldsymbol J-\boldsymbol I)\boldsymbol D\widehat{\boldsymbol p}},
#' where \eqn{\widehat{\boldsymbol p}=(\widehat{p}_0, \widehat{p}_1, \ldots, \widehat{p}_n)^\top},
#' with \eqn{\widehat{p}_k} being the standard estimate of the marginal probability for state
#' \eqn{s_k}, \eqn{\boldsymbol I} and \eqn{\boldsymbol J} are the identity and counteridentity
#' matrices of order \eqn{n + 1}, respectively, and \eqn{\boldsymbol D} is a pairwise distance
#' matrix for the elements in the set \eqn{\mathcal{S}} considering a specific distance
#' between ordinal states, \eqn{d(\cdot, \cdot)}. If \code{normalize = TRUE}, then the normalized estimate is computed, namely
#' \eqn{\frac{\widehat{asym}_{d}}{max_{s_i, s_j \in \mathcal{S}}d(s_i, s_j)}}.
#' @encoding UTF-8
#' @author
#' Ángel López-Oriona, José A. Vilar
#' @references{
#'
#'   \insertRef{weiss2019distance}{otsfeatures}
#'
#' }
#' @export

ordinal_asymmetry <- function(series, states, distance = 'Block',
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
  identity_matrix <- diag(n_states)
  counteridentity_matrix <- diag(rep(1, n_states))[n_states:1, ]
  d_0_m <- distance_function(states[1], states[n_states])

  distance_matrix <- matrix(0, nrow = n_states, ncol = n_states)

  for (i in 1 : n_states) {

    for (j in 1 : n_states) {

      distance_matrix[i, j] <- distance_function(states[i], states[j])

    }

  }

  estimated_asymmetry <- vector_mp %*% (counteridentity_matrix - identity_matrix) %*% distance_matrix %*% vector_mp
  normalized_estimated_asymmetry <- estimated_asymmetry/d_0_m

  if (normalize == FALSE ){

    return(estimated_asymmetry)

  } else {

    return(normalized_estimated_asymmetry)

  }

}
