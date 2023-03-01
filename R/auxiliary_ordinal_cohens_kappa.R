

auxiliary_ordinal_cohens_kappa_do1 <- function(series, states, max_lag = 10, alpha = 0.05) {

  series_length <- length(series)
  vector_cp <- c_marginal_probabilities(series, states)
  n_states <- length(states)
  coeff <- (series_length - 1)/series_length
  values_ordinal_cohens_kappa <- numeric(max_lag)
  estimated_dispersion <- coeff * ordinal_dispersion_2(series, states, distance = 'Block')

  for (i in 1 : max_lag) {

    values_ordinal_cohens_kappa[i] <- ordinal_cohens_kappa(series = series,
                                                           states = states,
                                                           distance = 'Block',
                                                           lag = i)

  }

  matrix_a_variance_2 <- matrix(0, nrow = n_states - 1, ncol = n_states - 1)

  for (i in 1 : (n_states - 1)) {

    for (j in 1 : (n_states - 1)) {

      matrix_a_variance_2[i, j] <- (vector_cp[min(i, j)] - vector_cp[i] * vector_cp[j])^2

    }

  }

  vector_test_statistic <- values_ordinal_cohens_kappa
  a_variance_1 <- 4/(series_length * estimated_dispersion^2)
  a_variance_2 <- sum(matrix_a_variance_2)
  a_variance <- a_variance_1 * a_variance_2
  a_sd <- sqrt(a_variance)
  vector_p_values <- 2 * (1 - stats::pnorm(abs(vector_test_statistic), mean = -1/series_length, sd = a_sd))
  critical_value <- stats::qnorm(1 - alpha/2, mean = -1/series_length, sd = a_sd)

  return_list <- list(values_ordinal_cohens_kappa = values_ordinal_cohens_kappa,
                      vector_p_values = vector_p_values,
                      critical_value = critical_value)

  return(return_list)

}
