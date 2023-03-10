#'
#' @title SyntheticData1
#' @description Synthetic dataset containing 80 OTS generated from four
#' different generating processes.
#' @usage data(SyntheticData1)
#' @format A \code{list} with two elements, which are:
#' \describe{
#' \item{\code{data}}{A list with 80 OTS.}
#' \item{\code{classes}}{A numeric vector indicating the corresponding classes associated with the elements in \code{data}.}
#' }
#' @details Each element in \code{data} is a 6-state OTS of length 600.
#' Series 1-20, 21-40, 41-60 and 61-80 were generated from
#' binomial AR(p) processes with different coefficients (see Scenario 1 in \insertCite{lopez2023submitted;textual}{otsfeatures}).
#' Therefore, there are 4 different classes in the dataset.
#' @references{
#'
#'   \insertRef{lopez2023submitted}{otsfeatures}
#'
#' }
"SyntheticData1"


