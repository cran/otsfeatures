#'
#' @title CreditRatings
#' @description Ordinal time series (OTS) of monthly credit ratings of different
#' European countries
#' @usage data(CreditRatings)
#' @format A \code{list} with one element, which is:
#' \describe{
#' \item{\code{data}}{A list with 28 MTS.}
#' }
#' @details Each element in \code{data} is an ordinal time series
#' containing 23 states (monthly credit ratings). The 28 countries of the European Union plus
#' the United Kingdom are considered. The sample period spans from January 2000 to December 2017, thus resulting serial realizations of length \eqn{T=216}.
#' For more information, see \insertCite{weiss2019distance;textual}{otsfeatures}.
#' @references{
#'
#'   \insertRef{weiss2019distance}{otsfeatures}
#'
#' }
"CreditRatings"


