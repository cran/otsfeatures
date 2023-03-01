

#' Computes the total mixed cumulative quantile correlation (TMCQC) between an ordinal and a
#' real-valued time series
#'
#' \code{total_mixed_c_correlation_2} returns the TMCQC
#' between an ordinal and a real-valued time series
#'
#' @param o_series An OTS.
#' @param n_series A real-valued time series.
#' @param lag The considered lag (default is 1).
#' @param states A numerical vector containing the corresponding
#' states.
#' @param features Logical. If \code{features = FALSE} (default), the value of the TMCLC is returned. Otherwise, the function
#' returns a vector with the individual components of the TMCQC.
#' @return If \code{features = FALSE} (default), returns the value of the TMCQC. Otherwise, the function
#' returns a vector of features, i.e., the vector contains the features employed to compute the
#' TMCLC.
#' @examples
#' tmclc <- total_mixed_c_correlation_2(o_series = SyntheticData1$data[[1]],
#' n_series = rnorm(600), states = 0 : 5) # Computing the TMCQC
#' # between the first series in dataset SyntheticData1 and white noise
#' feature_vector <- total_mixed_c_correlation_2(o_series = SyntheticData1$data[[1]],
#' n_series = rnorm(600), states = 0 : 5, features = TRUE) # Computing the corresponding
#' # vector of features
#' @details
#' Given a OTS of length \eqn{T} with range \eqn{\mathcal{S}=\{s_0, s_1, \ldots, s_n\}},
#' \eqn{\overline{X}_t=\{\overline{X}_1,\ldots, \overline{X}_T\}}, and
#' the cumulative binarized time series, which is defined as
#' \eqn{\overline{\boldsymbol Y}_t=\{\overline{\boldsymbol Y}_1, \ldots, \overline{\boldsymbol Y}_T\}},
#' with \eqn{\overline{\boldsymbol Y}_k=(\overline{Y}_{k,0}, \ldots, \overline{Y}_{k,n-1})^\top}
#' such that \eqn{\overline{Y}_{k,i}=1} if \eqn{\overline{X}_k \leq s_i} (\eqn{k=1,\ldots,T
#' , i=0,\ldots,n-1}), the function computes the estimated TMCQC given by
#' \deqn{\widehat{\Psi}_2^m(l)=\frac{1}{n}\sum_{i=0}^{n-1}\int_{0}^{1}\widehat{\psi}^\rho_{i}(l)^2d\rho,} where
#' \eqn{\widehat{\psi}_{i}^\rho(l)=\widehat{Corr}\big(Y_{t,i}, I(Z_{t-l}\leq q_{Z_t}(\rho)) \big)}, with
#' \eqn{\overline{Z}_t=\{\overline{Z}_1,\ldots, \overline{Z}_T\}} being a
#' \eqn{T}-length real-valued time series, \eqn{\rho \in (0, 1)} a probability
#' level, \eqn{I(\cdot)} the indicator function and \eqn{q_{Z_t}} the quantile
#' function of the corresponding real-valued process. If \code{features = TRUE}, the function
#' returns a vector whose components are the quantities \eqn{\int_{0}^{1}\widehat{\psi}^\rho_{i}(l)^2d\rho},
#' \eqn{i=0,1, \ldots,n-1}.
#' @encoding UTF-8
#' @author
#' Ángel López-Oriona, José A. Vilar
#' @export

total_mixed_c_correlation_2 <- function(o_series, n_series, lag = 1,
                                      states, features = FALSE) {
  check_ots(o_series)
  check_ots(n_series)
  series_length <- length(o_series)
  n_states <- length(states)
  grid_quantile <- seq(0, 1, by = 0.01)
  l_grid <- base::length(grid_quantile)
  c_binarized_series <- c_binarization(series = o_series, states = states)
  c_binarized_series_1 <- c_binarized_series[(lag + 1) : series_length,]
  c_binarized_series_2 <- c_binarized_series[1 : (series_length - lag),]

  correlation_matrix <- matrix(0, nrow = n_states - 1, ncol = l_grid)

  for (i in 1 : (n_states - 1)) {

    for (j in 1 : l_grid) {

      if (lag >= 0) {

        n_series_quantile <- as.numeric(n_series[1 : (series_length - lag)] <= grid_quantile[j])
        correlation_matrix[i, j] <- stats::cor(c_binarized_series_1[,i], n_series_quantile)^2

      } else {

        n_series_quantile <- as.numeric(n_series[(lag + 1) : series_length] <= grid_quantile[j])
        correlation_matrix[i, j] <- stats::cor(c_binarized_series_2[,i], n_series_quantile)^2

      }

    }

  }

  correlation_matrix[is.na(correlation_matrix)] <- 0


  vector_integrals <- numeric()

  for (i in 1 : (n_states - 1)) {

    vector_integrals[i] <- Bolstad2::sintegral(grid_quantile, correlation_matrix[i,])$int

  }


  if (features == FALSE) {

    return(mean(vector_integrals))

  } else {

    return(vector_integrals)

  }



}
