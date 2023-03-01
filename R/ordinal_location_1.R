

#' Computes the standard estimated location of an ordinal time series
#'
#' \code{ordinal_location_1} computes the standard estimated location
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
#' @param normalize Logical. If \code{normalize = FALSE} (default), the value of the standard estimated location is returned. Otherwise, the function
#' returns the normalized standard estimated location.
#' @return The standard estimated location.
#' @examples
#' estimated_location <- ordinal_location_1(series = AustrianWages$data[[100]],
#' states = 0 : 5) # Computing the standard location estimate
#' # for one series in dataset AustrianWages using the block distance
#' @details
#' Given an OTS of length \eqn{T} with range \eqn{\mathcal{S}=\{s_0, s_1, s_2, \ldots, s_n\}} (\eqn{s_0 < s_1 < s_2 < \ldots < s_n}),
#' \eqn{\overline{X}_t=\{\overline{X}_1,\ldots, \overline{X}_T\}}, the function computes the standard
#' estimated location given by \eqn{\widehat{x}_{loc, d}=}argmin\eqn{_{s \in \mathcal{S}}\frac{1}{T}\sum_{t=1}^Td\big(\overline{X}_t, s\big)},
#' where \eqn{d(\cdot, \cdot)} is a distance between ordinal states.
#' @encoding UTF-8
#' @author
#' Ángel López-Oriona, José A. Vilar
#' @references{
#'
#'   \insertRef{weiss2019distance}{otsfeatures}
#'
#' }
#' @export

ordinal_location_1 <- function(series, states, distance = 'Block',
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



    location_vector <- numeric(n_states)

    for (i in 1 : n_states) {

      distance_vector <- unlist(lapply(series,
                                       function(x) {distance_function(x, states[i])}))
      location_vector[i] <- base::mean(distance_vector)

    }

    min_position <- which.min(location_vector)

    if (normalize == FALSE ){

    return(states[min_position])

    } else {

    return(states[min_position]/(n_states - 1))

    }





}
