

#' Performs the hypothesis test associated with the
#' ordinal asymmetry for the block distance
#'
#' \code{test_ordinal_asymmetry} performs the hypothesis test associated with the
#' ordinal asymmetry for the block distance
#'
#' @param series An OTS (numerical vector with integers).
#' @param states A numeric vector containing the corresponding
#' states.
#' @param true_asymmetry The value for the true asymmetry.
#' @param alpha The significance level (default is 0.05).
#' @param temporal Logical. If \code{temporal = TRUE} (default), the test is performed for a time series. Otherwise,
#' the test is performed for i.i.d. data.
#' @param max_lag If \code{temporal = TRUE}, the maximum considered lag to compute the
#' estimates related to the cumulative joint probabilities.
#' @return The results of the hypothesis test.
#' @examples
#' results_test <- test_ordinal_asymmetry(AustrianWages$data[[100]],
#' states = 0 : 5, true_asymmetry = 2) # Performing the hypothesis test associated with the
#' # ordinal asymmetry for one OTS in dataset AustrianWages
#' @details
#' If \code{temporal = TRUE} (default), the function performs the hypothesis test based on the
#' ordinal asymmetry relying on Theorem 7.1.1 in \insertCite{weiss2019distance;textual}{otsfeatures}. Otherwise,
#' the test based on Theorem 4.1 in \insertCite{weiss2019distance;textual}{otsfeatures} is carried out.
#' @encoding UTF-8
#' @author
#' Ángel López-Oriona, José A. Vilar
#' @references{
#'
#'   \insertRef{weiss2019distance}{otsfeatures}
#'
#' }
#' @export

test_ordinal_asymmetry <- function(series, states, true_asymmetry,
                                    alpha = 0.05, temporal = TRUE, max_lag = 1) {

  check_ots(series)
  series_length <- length(series)
  n_states <- length(states)
  n <- n_states - 1
  test_statistic <- ordinal_asymmetry(series, states, distance = 'Block')
  estimated_dispersion <- ordinal_dispersion_2(series, states, distance = 'Block')
  vector_cp <- c_marginal_probabilities(series, states)
  coeff_1 <- 2/series_length
  coeff_2 <- 16/series_length

  vector_1 <- numeric()

  for (i in 1 : (n_states - 1)) {

      vector_1[i] <- (vector_cp[min(i, n_states - i)] - vector_cp[i] * vector_cp[n_states - i])

  }

  a_mean <- true_asymmetry + (1/series_length) * estimated_dispersion +
    coeff_1 * sum(vector_1)

  matrix_2 <- matrix(0, nrow = n_states - 1, ncol = n_states - 1)

  for (i in 1 : (n_states - 1)) {

    for (j in 1 : (n_states - 1)) {

      factor_1 <- 1 - vector_cp[i] - vector_cp[n_states - i]
      factor_2 <- 1 - vector_cp[j] - vector_cp[n_states - j]
      factor_3 <- vector_cp[min(i,j)] - vector_cp[i] * vector_cp[j]
      matrix_2[i, j] <- factor_1 * factor_2 * factor_3

    }

  }

  a_variance <- coeff_2 * sum(matrix_2)

  if (temporal == FALSE) {

    a_sd <- sqrt(a_variance)
    new_test_statistic <- (test_statistic - a_mean)/a_sd
    p_value <- 2 * (1 - stats::pnorm(abs(new_test_statistic)))
    critical_value <- stats::qnorm(1 - alpha/2)



    return_list <- list(test_statistic = new_test_statistic,
                        p_value = p_value,
                        critical_value = critical_value)

    return(return_list)

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

  a_mean_temporal <- a_mean + coeff_extra_1 * sum(vector_mean)

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

  a_variance_temporal <- a_variance + coeff_extra_2 * sum(vector_variance)

  if (temporal == TRUE) {

    a_sd_temporal <- sqrt(a_variance_temporal)
    new_test_statistic <- (test_statistic - a_mean_temporal)/a_sd_temporal
    p_value <- 2 * (1 - stats::pnorm(abs(new_test_statistic)))
    critical_value <- stats::qnorm(1 - alpha/2)



    return_list <- list(test_statistic = new_test_statistic,
                        p_value = p_value,
                        critical_value = critical_value)

    return(return_list)

  }





}
