

#' Computes the standard estimated dispersion of an ordinal time series
#'
#' \code{ordinal_dispersion_1} computes the standard estimated dispersion
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
#' @param normalize Logical. If \code{normalize = FALSE} (default), the value of the standard estimated dispersion is returned. Otherwise, the function
#' returns the normalized standard estimated dispersion.
#' @return The standard estimated dispersion.
#' @examples
#' estimated_dispersion <- ordinal_dispersion_1(series = AustrianWages$data[[100]],
#' states = 0 : 5) # Computing the standard dispersion estimate
#' # for one series in dataset AustrianWages using the block distance
#' @details
#' Given an OTS of length \eqn{T} with range \eqn{\mathcal{S}=\{s_0, s_1, s_2, \ldots, s_n\}} (\eqn{s_0 < s_1 < s_2 < \ldots < s_n}),
#' \eqn{\overline{X}_t=\{\overline{X}_1,\ldots, \overline{X}_T\}}, the function computes the standard
#' estimated dispersion given by \eqn{\widehat{disp}_{loc, d}=\frac{1}{T}\sum_{t=1}^Td\big(\overline{X}_t, \widehat{x}_{loc, d}\big)},
#' where \eqn{\widehat{x}_{loc, d}} is the standard estimate of the location and \eqn{d(\cdot, \cdot)} is a distance between ordinal states.
#' If \code{normalize = TRUE}, then the normalized dispersion is computed, namely
#' \eqn{\widehat{disp}_{loc, d}/}max\eqn{_{s_i, s_j \in \mathcal{S}}d(s_i, s_j)}.
#' @encoding UTF-8
#' @author
#' Ángel López-Oriona, José A. Vilar
#' @references{
#'
#'   \insertRef{weiss2019distance}{otsfeatures}
#'
#' }
#' @export

ordinal_dispersion_1 <- function(series, states, distance = 'Block',
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

  estimated_location <- ordinal_location_1(series, states, distance)
  distance_vector <- unlist(lapply(series,
                                   function(x) {distance_function(x, estimated_location)}))

  if (normalize == FALSE) {

    return(mean(distance_vector))

  } else {

    distance_matrix <- matrix(0, nrow = n_states, ncol = n_states)

    for (i in 1 : n_states) {


      for (j in (i + 1) : n_states) {

        if (i < n_states) {

        distance_matrix[i, j] <- distance_function(states[i], states[j])

        }

      }

    }

    return(mean(distance_vector)/max(distance_matrix))

  }



}
