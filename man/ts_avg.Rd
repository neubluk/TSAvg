% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ts_avg.R
\name{ts_avg}
\alias{ts_avg}
\title{ts_avg}
\usage{
ts_avg(
  x,
  model,
  k,
  h = 1,
  type = c("simple_avg", "simple_avg_n", "dist_avg_n", "dist_avg", "global_avg",
    "perf_avg", "perf_avg_r"),
  xtest_idx = NULL,
  benchmark = c("rw", "model"),
  ...
)
}
\arguments{
\item{x}{list of matrices. each matrix corresponds to one time series whereby the first column always corresponds to its index.}

\item{model}{a generic model function. the resulting model object must have forecast, refit, fitted, and residuals functions (as demonstrated in my_ets).}

\item{k}{vector of number of neighbors used for each time series in x.}

\item{h}{forecast horizon. if h>1, then multiple forecasts are computed. default=1.}

\item{type}{averaging type. one of  "simple_avg", "simple_avg_n", "dist_avg_n", "dist_avg", "global_avg", "perf_avg", "perf_avg_r".}

\item{xtest_idx}{index to split the every time series in train and test set. default = NULL.}

\item{benchmark}{benchmark forecast type for computing the errors. either random walk forecast "rw" or model forecast "model".}

\item{...}{additional parameters for model function}
}
\value{
object of class ts_avg
}
\description{
Perform TS Avg Forecasting for a particular type
}
\examples{
source("inst/my_ets.R")

x <- replicate(5,{
  len <- rpois(1,20)
  cbind(seq(rpois(1,5),length.out=len),matrix(rnorm(2*len),ncol = 2))
  },simplify=FALSE)


ts_avg_res <- ts_avg(x,my_ets,h=2,k=2,type="simple_avg", xtest_idx = 20)
plot(ts_avg_res)

ts_avg_res_2 <- ts_avg(lapply(x,"[",TRUE,1:2),my_ets,h=2,k=1,type="simple_avg", xtest_idx = 15)
plot(ts_avg_res_2,h=1)

# cv examples
res_mv <- ts_avg.cv(x,my_ets,min_ts=10,k_grid=2:5,xtest_idx = 20,verbose=TRUE)

plot(res_mv, "cv")
plot(res_mv, "error")

res_uv <- ts_avg.cv(lapply(x,"[",TRUE,1:2),my_ets,min_ts=10,k_grid=2:5, xtest_idx=20)

plot(res_uv, "best")
plot(res_uv, "cv")
}
\seealso{
\link{plot.ts_avg}
}
