

#' Performs the hypothesis test associated with the
#' ordinal dispersion for the block distance
#'
#' \code{test_ordinal_dispersion} performs the hypothesis test associated with the
#' ordinal dispersion for the block distance
#'
#' @param series An OTS (numerical vector with integers).
#' @param states A numeric vector containing the corresponding
#' states.
#' @param true_dispersion The value for the true dispersion.
#' @param alpha The significance level (default is 0.05).
#' @param temporal Logical. If \code{temporal = TRUE} (default), the test is performed for a time series. Otherwise,
#' the test is performed for i.i.d. data.
#' @param max_lag If \code{temporal = TRUE}, the maximum considered lag to compute the
#' estimates related to the cumulative joint probabilities.
#' @return The results of the hypothesis test.
#' @examples
#' results_test <- test_ordinal_dispersion(AustrianWages$data[[100]],
#' states = 0 : 5, true_dispersion = 2) # Performing the hypothesis test associated with the
#' # ordinal dispersion for one OTS in dataset AustrianWages
#' @details
#' If \code{temporal = TRUE} (default), the function performs the hypothesis test based on the
#' ordinal dispersion relying on Theorem 7.1.1 in \insertCite{weiss2019distance;textual}{otsfeatures}. Otherwise,
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

test_ordinal_dispersion <- function(series, states, true_dispersion,
                                  alpha = 0.05, temporal = TRUE, max_lag = 1) {

  check_ots(series)
  series_length <- length(series)
  n_states <- length(states)
  test_statistic <- ordinal_dispersion_2(series, states, distance = 'Block')
  a_mean <- (1 - 1/series_length) * true_dispersion
  vector_cp <- c_marginal_probabilities(series, states)
  coeff <- 4/series_length

  matrix_1 <- matrix(0, nrow = n_states - 1, ncol = n_states - 1)

  for (i in 1 : (n_states - 1)) {

    for (j in 1 : (n_states - 1)) {

      matrix_1[i, j] <- (1 - 2 * vector_cp[i]) * (1 - 2 * vector_cp[j]) * (vector_cp[min(i,j)] - vector_cp[i] * vector_cp[j])

    }

  }

  a_variance <- coeff * sum(matrix_1)

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

  coeff_extra <- 8/series_length

  vector_mean <- numeric()

  for (i in 1 : max_lag) {

    factor_2_prev_1 <- diag(c_joint_probabilities(series, i, states))
    factor_2_prev_2 <- vector_cp^2
    factor_2 <- sum(factor_2_prev_1 - factor_2_prev_2)
    vector_mean[i] <- factor_2

  }

  a_mean_temporal <- a_mean - coeff * sum(vector_mean)

  vector_variance <- numeric()

  for (i in 1 : max_lag) {

    factor_2_prev <- matrix(0, nrow = n_states - 1,
                              ncol = n_states - 1)

    for (j in 1 : (n_states - 1)) {

      for (k in 1 : (n_states - 1)) {

        matrix_jp <- c_joint_probabilities(series, i, states)
        factor_2_prev[j, k] <- (1 - 2 * vector_cp[j]) * (1 - 2 * vector_cp[k]) * (matrix_jp[j,k] - vector_cp[j] * vector_cp[k])

      }

    }

    vector_variance[i] <- sum(factor_2_prev)

  }

  a_variance_temporal <- a_variance + coeff_extra * sum(vector_variance)

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
