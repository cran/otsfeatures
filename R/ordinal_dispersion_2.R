

#' Computes the estimated dispersion of an ordinal time series according to
#' the approach based on the diversity coefficient (DIVC)
#'
#' \code{ordinal_dispersion_2} computes the estimated dispersion
#' of an ordinal time series according to the approach based on the
#' diversity coefficient
#'
#' @param series An OTS.
#' @param states A numerical vector containing the corresponding
#' states.
#' @param distance A function defining the underlying distance between
#' states. The Hamming, block and Euclidean distances are already
#' implemented by means of the arguments "Hamming", "Block" (default)
#' and "Euclidean". Otherwise, a function taking as input two states must
#' be provided.
#' @param normalize Logical. If \code{normalize = FALSE} (default), the value of the estimated dispersion is returned. Otherwise, the function
#' returns the normalized estimated dispersion.
#' @return The estimated dispersion according to the approach based on the
#' diversity coefficient.
#' @examples
#' estimated_dispersion <- ordinal_dispersion_2(series = AustrianWages$data[[100]],
#' states = 0 : 5) # Computing the DIVC dispersion estimate
#' # for one series in dataset AustrianWages using the block distance
#' @details
#' Given an OTS of length \eqn{T} with range \eqn{\mathcal{S}=\{s_0, s_1, s_2, \ldots, s_n\}} (\eqn{s_0 < s_1 < s_2 < \ldots < s_n}),
#' \eqn{\overline{X}_t=\{\overline{X}_1,\ldots, \overline{X}_T\}}, the function computes the DIVC
#' estimated dispersion given by \eqn{\widehat{disp}_{d}=\frac{T}{T-1}\sum_{i,j=0}^nd\big(s_i, s_j\big)\widehat{p}_i\widehat{p}_j},
#' where \eqn{d(\cdot, \cdot)} is a distance between ordinal states and \eqn{\widehat{p}_k} is the
#' standard estimate of the marginal probability for state \eqn{s_k}.
#' If \code{normalize = TRUE}, and \code{distance = "Block"} or \code{distance = "Euclidean"}, then the normalized versions are computed, that is,
#' the corresponding estimates are divided by the factors \eqn{2/m} or \eqn{2/m^2}, respectively.
#' @encoding UTF-8
#' @author
#' Ángel López-Oriona, José A. Vilar
#' @references{
#'
#'   \insertRef{weiss2019distance}{otsfeatures}
#'
#' }
#' @export

ordinal_dispersion_2 <- function(series, states, distance = 'Block',
                                 normalize = FALSE) {

  check_ots(series)
  series_length <- length(series)
  coeff <- series_length/(series_length - 1)
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
  distance_matrix <- matrix(0, nrow = n_states, ncol = n_states)

  for (i in 1 : n_states) {

    for (j in 1 : n_states) {

        distance_matrix[i, j] <- distance_function(states[i], states[j])

    }

  }

  nonnormalized_dispersion <- as.numeric(coeff * vector_mp %*% distance_matrix %*% vector_mp)

  if (normalize == FALSE) {

    return(nonnormalized_dispersion)

  } else {

    if (distance == 'Block' | distance == 'Euclidean') {

      if (distance == 'Block') {

        return(4/(n_states - 1)*nonnormalized_dispersion)

      } else {

        return(2/(n_states - 1)^2*nonnormalized_dispersion)

      }

    } else {

      return(nonnormalized_dispersion)

    }

  }



}
