% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{cif2016}
\alias{cif2016}
\title{CIF 2016 Competition Dataset}
\format{
A list containing
\enumerate{
\item a long data.frame with series_name, horizon, and series_value
\item frequency
\item forecast horizon if not given in the data.frame
\item contains missing values
\item contains time series of equal length
}
}
\source{
\url{https://zenodo.org/record/4656042}
}
\usage{
cif2016
}
\description{
The Computational Intelligence in Forecasting (CIF) 2016 dataset contains 72 monthly time series originated from the banking domain used in the CIF 2016 forecasting competition. Out of the 72 time series, 24 series contains real-time data while the remaining 48 series are artificially generated. There are two forecast horizons considered in the competition where 57 series consider the forecast horizon as 12 and the remaining 15 series consider the forecast horizon as 6.
}
\references{
Stepnicka, M., Burda, M., 2017. On the results and observations of the time series forecasting competition CIF 2016. In: 2017 IEEE International Conference on Fuzzy Systems (FUZZ-IEEE). pp. 1–6.
}
\keyword{datasets}
