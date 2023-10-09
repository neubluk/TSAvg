# TSAvg

[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

This R package features the essential functions from [Neubauer and Filzmoser (2023)](<https://arxiv.org/abs/2306.07119>). Exported functions are

-   `ts_avg`: performs TS Averaging for forecasting
-   `ts_avg.cv`: performs TSCV for selecting the optimal number of neighbors per time series, and yields the best models.

## Installation

``` s
# install.packages("remotes")
remotes::install_github("neubluk/TSAvg")
```

## Usage

The functionality requires a model as in `my_ets.R` which yields an object with `forecast`, `refit`, `fitted`, and `residuals` functions.

``` s
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
plot(res_mv, "best")
```

## License

This package is free and open source software, licensed under GPL-3.
