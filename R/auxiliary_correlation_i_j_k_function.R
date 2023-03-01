

correlation_i_j_k_function <- function(series, i_stat, j_stat, k = 1, states) {

  series_length <- length(series)
  c_binarized_series <- c_binarization(series, states)
  a <- c_binarized_series[(k + 1) : series_length, i_stat + 1]
  b <- c_binarized_series[1 : (series_length - k), j_stat + 1]


  return(stats::cor(a, b))

}
