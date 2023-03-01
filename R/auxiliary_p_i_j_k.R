

p_i_j_k_function <- function(series, i_stat, j_stat, k = 1) {

  series_length <- length(series)
  a <- series[(k + 1) : series_length]
  b <- series[1 : (series_length - k)]


  number <- series_length - k
  count <- numeric(number)

  for (i in 1 : number) {

    if (a[i] == i_stat & b[i] == j_stat) {

      count[i] <- 1

    } else {

      count[i] <- 0

    }

  }


  return(sum(count)/(series_length-k))

}
