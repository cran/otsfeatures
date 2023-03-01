

#' Constructs a confidence interval for the
#' ordinal dispersion (block distance)
#'
#' \code{ci_ordinal_dispersion} constructs a confidence interval for the
#' ordinal dispersion (block distance)
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
#' ci_dispersion <- ci_ordinal_dispersion(AustrianWages$data[[100]],
#' states = 0 : 5) # Constructing a confidence interval for the
#' # ordinal dispersion for one OTS in dataset AustrianWages
#' @details
#' If \code{temporal = TRUE} (default), the function constructs the confidence interval for the
#' ordinal dispersion relying on Theorem 7.1.1 in \insertCite{weiss2019distance;textual}{otsfeatures}. Otherwise,
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

ci_ordinal_dispersion <- function(series, states, level = 0.95, temporal = TRUE, max_lag = 1) {

  check_ots(series)
  alpha <- 1 - level
  series_length <- length(series)
  n_states <- length(states)
  estimated_dispersion <- ordinal_dispersion_2(series, states, distance = 'Block')
  vector_cp <- c_marginal_probabilities(series, states)
  coeff <- 4/series_length
  quantile_normal <- stats::qnorm(1 - alpha/2)
  denominator <- 1 - 1/series_length

  matrix_1 <- matrix(0, nrow = n_states - 1, ncol = n_states - 1)

  for (i in 1 : (n_states - 1)) {

    for (j in 1 : (n_states - 1)) {

      matrix_1[i, j] <- (1 - 2 * vector_cp[i]) * (1 - 2 * vector_cp[j]) * (vector_cp[min(i,j)] - vector_cp[i] * vector_cp[j])

    }

  }

  a_variance <- coeff * sum(matrix_1)
  a_sd <- sqrt(a_variance)

  if (temporal == FALSE) {

    lower_bound <- (estimated_dispersion - quantile_normal * a_sd)/denominator
    upper_bound <- (estimated_dispersion + quantile_normal * a_sd)/denominator
    return_df <- data.frame(lower_bound, upper_bound)
    colnames(return_df) <- c('Lower bound', 'Upper bound')
    return(return_df)

  }

  coeff_extra <- 8/series_length

  vector_mean <- numeric()

  for (i in 1 : max_lag) {

    factor_2_prev_1 <- diag(c_joint_probabilities(series, i, states))
    factor_2_prev_2 <- vector_cp^2
    factor_2 <- sum(factor_2_prev_1 - factor_2_prev_2)
    vector_mean[i] <- factor_2

  }

  coeff_mean <- coeff * sum(vector_mean)

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
  a_sd_temporal <- sqrt(a_variance_temporal)

  if (temporal == TRUE) {

    lower_bound <- (estimated_dispersion - quantile_normal * a_sd_temporal + coeff_mean)/denominator
    upper_bound <- (estimated_dispersion + quantile_normal * a_sd_temporal + coeff_mean)/denominator
    return_df <- data.frame(lower_bound, upper_bound)
    colnames(return_df) <- c('Lower bound', 'Upper bound')
    return(return_df)

  }



}
