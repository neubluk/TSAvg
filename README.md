# TSAvg

[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

This R package features the essential functions from [Neubauer and Filzmoser (2023)](<https://doi.org/10.1016/j.ijforecast.2024.02.002>). Exported functions are

-   `ts_avg`: performs TS Averaging for forecasting
-   `ts_avg.cv`: performs TSCV for selecting the optimal number of neighbors per time series, and yields the best models.

## Citation

```s
@article{NEUBAUER2024,
  title = {Improving forecasts for heterogeneous time series by “averaging”, with application to food demand forecasts},
  journal = {International Journal of Forecasting},
  year = {2024},
  issn = {0169-2070},
  doi = {https://doi.org/10.1016/j.ijforecast.2024.02.002},
  url = {https://www.sciencedirect.com/science/article/pii/S0169207024000074},
  author = {Lukas Neubauer and Peter Filzmoser},
  keywords = {Time series, Forecasting, Combining forecasts, Dynamic time warping, k-nearest neighbors},
  abstract = {A common forecasting setting in real-world applications considers a set of possibly heterogeneous time series of the same domain. Due to the different properties of each time series, such as length, obtaining forecasts for each individual time series in a straightforward way is challenging. This paper proposes a general framework utilizing a similarity measure in dynamic time warping to find similar time series to build neighborhoods in a k-nearest neighbor fashion and improve forecasts of possibly simple models by averaging. Several ways of performing the averaging are suggested, and theoretical arguments underline the usefulness of averaging for forecasting. Additionally, diagnostic tools are proposed for a deep understanding of the procedure.}
}
```

## Installation

``` s
# install.packages("remotes")
remotes::install_github("neubluk/TSAvg")
```

## Usage

The functionality requires a model as in `demo/my_ets.R` which yields an object with `forecast`, `refit`, `fitted`, and `residuals` functions. A list containing all of those functions is needed as an input for the *TSAvg* functions.

``` s
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
```

## Data Examples

Data examples are available in `demo/data_examples.R`.

## Data

Two datasets are available.

- Food Demand Data from Schrankel GmbH (```food_demand```)
- CIF 2016 Data, taken from Stepnicka, M., Burda, M., 2017. On the results and observations of the time series forecasting competition CIF 2016. In: 2017 IEEE International Conference on Fuzzy Systems (FUZZ-IEEE). pp. 1–6. (```cif2016```)

## License

This package is free and open source software, licensed under GPL-3.
