% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ts_avg.R
\name{plot.ts_avg}
\alias{plot.ts_avg}
\title{plot.ts_avg}
\usage{
\method{plot}{ts_avg}(
  object,
  plot_type = c("response", "error", "error_quants"),
  h = 1,
  ids = NULL,
  ref_plot = NULL
)
}
\arguments{
\item{object}{result of ts_avg function}

\item{plot_type}{response or error}

\item{h}{forecast horizon to use. default=1.}

\item{ids}{vector of ids to filter appropriate plots.}

\item{ref_plot}{ts_avg plots can be concatenated such that different types and ks can be shwon together.}
}
\value{
ggplot2 object
}
\description{
Plot ts_avg result
}
\examples{
source("demo/my_ets.R")

x <- replicate(100,{
  len <- rpois(1,20)
  cbind(seq(rpois(1,5),length.out=len),matrix(rnorm(2*len),ncol = 2))
  },simplify=FALSE)

model <- list(model = my_ets,
              forecast = forecast.my_ets,
              refit = refit.my_ets,
              fitted = fitted.my_ets,
              residuals = residuals.my_ets)

ts_avg_res <- ts_avg(
  x,
  model,
  h = 4,
  k = 2,
  type = "global_avg",
  xtest_idx = 20,
  execution_time = TRUE,
  detailed_result = FALSE
)
plot(ts_avg_res)

ts_avg_res_2 <- ts_avg(lapply(x,"[",TRUE,1:2),model,h=2,k=1,type="simple_avg", xtest_idx = 15)
plot(ts_avg_res_2,h=1)

# cv examples
min_ts <- max(unlist(lapply(x, function(xi)min(xi[,1])+3)))
res_mv <- ts_avg.cv(
  x,
  model,
  min_ts = min_ts + 1,
  k_grid = c(1, 10, 20),
  xtest_idx = 20,
  verbose = TRUE
)

plot(res_mv, "cv")
plot(res_mv, "error")
plot(res_mv, "best")

p <- plot(res_mv$final_model, "error_quants")
plot(ts_avg_res, "error_quants", ref_plot=p)
}
