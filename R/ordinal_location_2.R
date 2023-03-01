

#' Computes the estimated location of an ordinal time series
#' with respect to the lowest category
#'
#' \code{ordinal_location_2} computes the estimated location
#' of an ordinal time series with respect to the lowest category
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
#' @return The estimated location with respect to the lowest category.
#' @examples
#' estimated_location <- ordinal_location_2(series = AustrianWages$data[[100]],
#' states = 0 : 5) # Computing the location estimate
#' # with respect to the lowest state for one series in dataset AustrianWages
#' @details
#' Given an OTS of length \eqn{T} with range \eqn{\mathcal{S}=\{s_0, s_1, s_2, \ldots, s_n\}} (\eqn{s_0 < s_1 < s_2 < \ldots < s_n}),
#' \eqn{\overline{X}_t=\{\overline{X}_1,\ldots, \overline{X}_T\}}, the function computes the
#' estimated location with respect to the lowest state, that is, the state
#' \eqn{s_j} such that \eqn{a_j=d(s_j, s_0)} is the closest to
#' \eqn{\frac{1}{T}\sum_{t=1}^Td\big(\overline{X}_t, s_0\big)} is determined,
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

ordinal_location_2 <- function(series, states, distance = 'Block',
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



    distance_vector <- unlist(lapply(series,
                                     function(x) {distance_function(x, states[1])}))
    mean_distance <- base::mean(distance_vector)
    vector_distances <- unlist(lapply(states,
                                      function(x) {distance_function(x, states[1])}))
    position <- which.min(abs(vector_distances - mean_distance))

    if (normalize == FALSE ){

      return(states[position])

    } else {

      return(states[position]/(n_states - 1))

    }





}
