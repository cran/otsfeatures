

#' Computes the total cumulative correlation of an ordinal time series
#'
#' \code{total_c_correlation} returns the value of the total cumulative correlation for
#' an ordinal time series
#'
#' @param series An OTS.
#' @param lag The considered lag (default is 1).
#' @param states A numerical vector containing the corresponding
#' states.
#' @param features Logical. If \code{features = FALSE} (default), the value of the total cumulative correlation is returned. Otherwise, the function
#' returns a matrix with the individual components of the total cumulative correlation
#' @return If \code{features = FALSE} (default), returns the value of the total cumulative correlation. Otherwise, the function
#' returns a matrix of features, i.e., the matrix contains the features employed to compute the
#' total cumulative correlation.
#' @examples
#' tcc <- total_c_correlation(series = AustrianWages$data[[100]],
#' states = 0 : 5) # Computing the total cumulative correlation
#' # for one of the series in dataset AustrianWages
#' feature_matrix <- total_c_correlation(series = AustrianWages$data[[100]],
#' states = 0 : 5) # Computing the corresponding matrix of features
#' @details
#' Given an OTS of length \eqn{T} with range \eqn{\mathcal{S}=\{s_0, s_1, \ldots, s_n\}},
#' \eqn{\overline{X}_t=\{\overline{X}_1,\ldots, \overline{X}_T\}}, and
#' the cumulative binarized time series, which is defined as
#' \eqn{\overline{\boldsymbol Y}_t=\{\overline{\boldsymbol Y}_1, \ldots, \overline{\boldsymbol Y}_T\}},
#' with \eqn{\overline{\boldsymbol Y}_k=(\overline{Y}_{k,0}, \ldots, \overline{Y}_{k,n-1})^\top}
#' such that \eqn{\overline{Y}_{k,i}=1} if \eqn{\overline{X}_k\leq s_i} (\eqn{k=1,\ldots,T,
#' , i=0,\ldots,n-1}), the function computes the estimated average \eqn{\widehat{\Psi}(l)^c=\frac{1}{n^2}\sum_{i,j=0}^{n-1}\widehat{\psi}_{ij}(l)^2},
#' where \eqn{\widehat{\psi}_{ij}(l)} is the estimated correlation
#' \eqn{\widehat{Corr}(Y_{t, i}, Y_{t-l, j})}, \eqn{i,j=0, 1,\ldots,n-1}. If \code{features = TRUE}, the function
#' returns a matrix whose components are the quantities \eqn{\widehat{\psi}_{ij}(l)},
#' \eqn{i,j=0,1, \ldots,n-1}.
#' @encoding UTF-8
#' @author
#' Ángel López-Oriona, José A. Vilar
#' @export


total_c_correlation <- function(series, lag = 1, states, features = FALSE) {

  check_ots(series)
  n_states <- length(states)
  matrix_prev <- base::matrix(0, nrow = n_states - 1, ncol = n_states - 1)

  for (i in 0 : (n_states - 2)) {

    for (j in 0 : (n_states - 2)) {

      matrix_prev[i + 1, j + 1] <- correlation_i_j_k_function(series = series,
                                                      i_stat = states[i + 1],
                                                      j_stat = states[j + 1],
                                                      k = lag,
                                                      states = states)

    }

  }

  matrix_prev[is.na(matrix_prev)] <- 0

  matrix_sq <- matrix_prev^2

  if (features == FALSE) {

    return(sum(matrix_sq))

  } else {

    return(matrix_prev)

  }

}
