

auxiliary_test_skewness <- function(series, states, lag = 1){

  vector_cp <- c_marginal_probabilities(series, states)
  matrix_cp_product <- vector_cp %*% t(vector_cp)
  matrix_cp_lag <- c_joint_probabilities(series, lag, states)

  return(sum(matrix_cp_lag - matrix_cp_product))

}
