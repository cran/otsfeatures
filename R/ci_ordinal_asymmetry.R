

#' Constructs a confidence interval for the
#' ordinal asymmetry (block distance)
#'
#' \code{ci_ordinal_asymmetry} constructs a confidence interval for the
#' ordinal asymmetry (block distance)
#'
#' @param series An OTS (numerical vector with integers).
#' @param states A numeric vector containing the corresponding
#' states.
#' @param level The confidence level (default is 0.95).
#' @param temporal Logical. If \code{temporal = TRUE} (default), the interval is computed for a time series. Otherwise,
#' the interval is computed for i.i.d. data.
#' @param max_lag If \code{temporal = TRUE}, the maximum considered lag to compute the
#' estimates related to the cumulative joint probabilities.
#' @return The confidence interval.
#' @examples
#' ci_asymmetry <- ci_ordinal_asymmetry(AustrianWages$data[[100]],
#' states = 0 : 5) # Constructing a confidence interval for the
#' # ordinal asymmetry for one OTS in dataset AustrianWages
#' @details
#' If \code{temporal = TRUE} (default), the function constructs the confidence interval for the
#' ordinal asymmetry relying on Theorem 7.1.1 in \insertCite{weiss2019distance;textual}{otsfeatures}. Otherwise,
#' the interval is constructed according to Theorem 4.1 in \insertCite{weiss2019distance;textual}{otsfeatures}.
#' @encoding UTF-8
#' @author
#' Ángel López-Oriona, José A. Vilar
#' @references{
#'
#'   \insertRef{weiss2019distance}{otsfeatures}
#'
#' }
#' @export

ci_ordinal_asymmetry <- function(series, states, level = 0.95, temporal = TRUE, max_lag = 1) {

  check_ots(series)
  alpha <- 1 - level
  series_length <- length(series)
  n_states <- length(states)
  estimated_dispersion <- ordinal_dispersion_2(series, states, distance = 'Block')
  estimated_asymmetry <- ordinal_asymmetry(series, states, distance = 'Block')
  vector_cp <- c_marginal_probabilities(series, states)
  coeff_1 <- 2/series_length
  coeff_2 <- 16/series_length
  quantile_normal <- stats::qnorm(1 - alpha/2)
  factor_a <- 1/series_length * estimated_dispersion

  vector_1 <- numeric()

  for (i in 1 : (n_states - 1)) {

    vector_1[i] <- (vector_cp[min(i, n_states - i)] - vector_cp[i] * vector_cp[n_states - i])

  }

  factor_b <- coeff_1 * sum(vector_1)

  matrix_2 <- matrix(0, nrow = n_states - 1, ncol = n_states - 1)

  for (i in 1 : (n_states - 1)) {

    for (j in 1 : (n_states - 1)) {

      factor_1 <- 1 - vector_cp[i] - vector_cp[n_states - i]
      factor_2 <- 1 - vector_cp[j] - vector_cp[n_states - j]
      factor_3 <- vector_cp[min(i,j)] - vector_cp[i] * vector_cp[j]
      matrix_2[i, j] <- factor_1 * factor_2 * factor_3

    }

  }

  factor_d <- coeff_2 * sum(matrix_2)

  if (temporal == FALSE) {

    lower_bound <- estimated_asymmetry - quantile_normal * factor_d - factor_a - factor_b
    upper_bound <- estimated_asymmetry + quantile_normal * factor_d - factor_a - factor_b
    return_df <- data.frame(lower_bound, upper_bound)
    colnames(return_df) <- c('Lower bound', 'Upper bound')
    return(return_df)

  }

  coeff_extra_1 <- 4/series_length
  coeff_extra_2 <- 32/series_length

  vector_mean <- numeric()

  for (i in 1 : max_lag) {

    factor_2_prev <- matrix(0, nrow = n_states - 1,
                            ncol = n_states - 1)

    for (j in 1 : (n_states - 1)) {

      for (k in 1 : (n_states - 1)) {

        matrix_jp <- c_joint_probabilities(series, i, states)
        factor_2_prev[j, k] <- matrix_jp[j, j] - vector_cp[j]^2 + matrix_jp[j, n_states - j] - vector_cp[j] * vector_cp[n_states - j]

      }

    }

    vector_mean[i] <- sum(factor_2_prev)

  }

  factor_c <- coeff_extra_1 * sum(vector_mean)


  vector_variance <- numeric()

  for (i in 1 : max_lag) {

    factor_2_prev <- matrix(0, nrow = n_states - 1,
                            ncol = n_states - 1)

    for (j in 1 : (n_states - 1)) {

      for (k in 1 : (n_states - 1)) {

        matrix_jp <- c_joint_probabilities(series, i, states)
        factor_2_prev[j, k] <- (1 - vector_cp[j] - vector_cp[n_states - j]) * (1 - vector_cp[k] -vector_cp[n_states - k]) * (matrix_jp[j, k] - vector_cp[j] * vector_cp[k])
      }

    }

    vector_variance[i] <- sum(factor_2_prev)

  }

factor_e <- coeff_extra_2 * sum(vector_variance)

if (temporal == TRUE) {

  lower_bound <- estimated_asymmetry - quantile_normal * sqrt(factor_d + factor_e) - factor_a - factor_b - factor_c
  upper_bound <- estimated_asymmetry + quantile_normal * sqrt(factor_d + factor_e) - factor_a - factor_b - factor_c
  return_df <- data.frame(lower_bound, upper_bound)
  colnames(return_df) <- c('Lower bound', 'Upper bound')
  return(return_df)

}



}
