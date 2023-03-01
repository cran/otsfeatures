

#' Computes the total mixed cumulative linear correlation (TMCLC) between an ordinal and a
#' real-valued time series
#'
#' \code{total_mixed_c_correlation_1} returns the TMCLC between an ordinal and a
#' real-valued time series
#'
#' @param o_series An OTS.
#' @param n_series A real-valued time series.
#' @param lag The considered lag (default is 1).
#' @param states A numerical vector containing the corresponding
#' states.
#' @param features Logical. If \code{features = FALSE} (default), the value of the TMCLC is returned. Otherwise, the function
#' returns a vector with the individual components of the TMCLC.
#' @return If \code{features = FALSE} (default), returns the value of the TMCLC. Otherwise, the function
#' returns a vector of features, i.e., the vector contains the features employed to compute the
#' TMCLC.
#' @examples
#' tmclc <- total_mixed_c_correlation_1(o_series = SyntheticData1$data[[1]],
#' n_series = rnorm(600), states = 0 : 5) # Computing the TMCLC
#' # between the first series in dataset SyntheticData1 and white noise
#' feature_vector <- total_mixed_c_correlation_1(o_series = SyntheticData1$data[[1]],
#' n_series = rnorm(600), states = 0 : 5, features = TRUE) # Computing the corresponding
#' # vector of features
#' @details
#' Given a OTS of length \eqn{T} with range \eqn{\mathcal{S}=\{s_0, s_1, \ldots, s_n\}},
#' \eqn{\overline{X}_t=\{\overline{X}_1,\ldots, \overline{X}_T\}}, and
#' the cumulative binarized time series, which is defined as
#' \eqn{\overline{\boldsymbol Y}_t=\{\overline{\boldsymbol Y}_1, \ldots, \overline{\boldsymbol Y}_T\}},
#' with \eqn{\overline{\boldsymbol Y}_k=(\overline{Y}_{k,0}, \ldots, \overline{Y}_{k,n-1})^\top}
#' such that \eqn{\overline{Y}_{k,i}=1} if \eqn{\overline{X}_k \leq s_i} (\eqn{k=1,\ldots,T
#' , i=0,\ldots,n-1}), the function computes the estimated TMCLC given by
#' \deqn{\widehat{\Psi}_1^m(l)=\frac{1}{n}\sum_{i=0}^{n-1}\widehat{\psi}_{i}^*(l)^2,} where
#' \eqn{\widehat{\psi}_{i}^*(l)=\widehat{Corr}(Y_{t,i}, Z_{t-l})}, with
#' \eqn{\overline{Z}_t=\{\overline{Z}_1,\ldots, \overline{Z}_T\}} being a
#' \eqn{T}-length real-valued time series. If \code{features = TRUE}, the function
#' returns a vector whose components are the quantities \eqn{\widehat{\psi}_{i}(l)},
#' \eqn{i=0,1, \ldots,n-1}.
#' @encoding UTF-8
#' @author
#' Ángel López-Oriona, José A. Vilar
#' @export

total_mixed_c_correlation_1 <- function(o_series, n_series, lag = 1,
                                      states, features = FALSE) {
  check_ots(o_series)
  check_ots(n_series)
  series_length <- length(o_series)
  n_states <- length(states)
  c_binarized_series <- c_binarization(series = o_series, states = states)
  c_binarized_series_1 <- c_binarized_series[(lag + 1) : series_length,]
  c_binarized_series_2 <- c_binarized_series[1 : (series_length - lag),]
  n_series_1 <- n_series[(lag + 1) : series_length]
  n_series_2 <- n_series[1 : (series_length - lag)]

  correlation_measure <- numeric(n_states - 1)

  for (i in 1 : (n_states - 1)) {

    if (lag >= 0) {

      correlation_measure[i] <- stats::cor(c_binarized_series_1[,i], n_series_2)

    } else {

      correlation_measure[i] <- stats::cor(c_binarized_series_2[,i], n_series_1)

    }

  }

  correlation_measure[is.na(correlation_measure)] <- 0


  sum_sq_correlation_measure <- sum(correlation_measure^2)

  if (features == FALSE) {

    return(sum_sq_correlation_measure/(n_states-1))

  } else {

    return(correlation_measure)

  }



}
