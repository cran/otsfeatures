

#' Computes the estimated ordinal Cohen's kappa of an ordinal time series
#'
#' \code{ordinal_cohens_kappa} computes the estimated ordinal Cohen's kappa
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
#' @param lag The considered lag.
#' @return The estimated ordinal Cohen's kappa.
#' @examples
#' estimated_ock <- ordinal_cohens_kappa(series = AustrianWages$data[[100]],
#' states = 0 : 5) # Computing the estimated ordinal Cohen's kappa
#' # for one series in dataset AustrianWages using the block distance
#' @details
#' Given an OTS of length \eqn{T} with range \eqn{\mathcal{S}=\{s_0, s_1, s_2, \ldots, s_n\}} (\eqn{s_0 < s_1 < s_2 < \ldots < s_n}),
#' \eqn{\overline{X}_t=\{\overline{X}_1,\ldots, \overline{X}_T\}}, the function computes the
#' estimated ordinal Cohen's kappa given by \eqn{\widehat{\kappa}_d(l)=\frac{\widehat{disp}_d(X_t)-\widehat{E}[d(X_t, X_{t-l})]}{{\widehat{disp}}_d(X_t)}},
#' where \eqn{\widehat{disp}_{d}(X_t)=\frac{T}{T-1}\sum_{i,j=0}^nd\big(s_i, s_j\big)\widehat{p}_i\widehat{p}_j} is the DIVC estimate of the dispersion, with
#' \eqn{d(\cdot, \cdot)} being a distance between ordinal states and \eqn{\widehat{p}_k} being the
#' standard estimate of the marginal probability for state \eqn{s_k},
#' and \eqn{\widehat{E}[d(X_t, X_{t-l})]=\frac{1}{T-l} \sum_{t=l+1}^T d(\overline{X}_t, \overline{X}_{t-l})}.
#' @encoding UTF-8
#' @author
#' Ángel López-Oriona, José A. Vilar
#' @references{
#'
#'   \insertRef{weiss2019distance}{otsfeatures}
#'
#' }
#' @export

ordinal_cohens_kappa <- function(series, states, distance = 'Block', lag = 1) {

  check_ots(series)
  series_length <- length(series)
  coeff <- (series_length - 1)/series_length
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

  estimated_dispersion <- coeff * ordinal_dispersion_2(series, states, distance = distance)
  vector_estimated_expectation <- numeric()

  for (i in (lag + 1) : series_length) {

    vector_estimated_expectation[i - (lag)] <- distance_function(series[i], series[i-lag])

  }

  estimated_expectation <- mean(vector_estimated_expectation)

  return((estimated_dispersion - estimated_expectation)/estimated_dispersion)



}
