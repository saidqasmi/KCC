% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estim_mar2.R
\name{Sigma_ar}
\alias{Sigma_ar}
\title{Compute the covariance matrix corresponding to an autoregressive processes of
order 1}
\usage{
Sigma_ar(alpha, n)
}
\arguments{
\item{alpha}{a numeric corresponding to the AR parameter at lag 1}

\item{n}{the whished number of lines/columns in the covariance matrix,
usually corresponds to the number of years of observations.}
}
\value{
a symmetric matrix corresponding to the AR process.
}
\description{
Compute the covariance matrix corresponding to an autoregressive processes of
order 1
}
\examples{

}
