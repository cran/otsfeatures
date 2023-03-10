% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SyntheticData2.R
\docType{data}
\name{SyntheticData2}
\alias{SyntheticData2}
\title{SyntheticData2}
\format{
A \code{list} with two elements, which are:
\describe{
\item{\code{data}}{A list with 80 OTS.}
\item{\code{classes}}{A numeric vector indicating the corresponding classes associated with the elements in \code{data}.}
}
}
\usage{
data(SyntheticData2)
}
\description{
Synthetic dataset containing 80 OTS generated from four
different generating processes.
}
\details{
Each element in \code{data} is a 6-state OTS of length 600.
Series 1-20, 21-40, 41-60 and 61-80 were generated from
binomial INARCH(p) processes with different coefficients (see Scenario 2 in \insertCite{lopez2023submitted;textual}{otsfeatures}).
Therefore, there are 4 different classes in the dataset.
}
\references{
{

  \insertRef{lopez2023submitted}{otsfeatures}

}
}
\keyword{datasets}
