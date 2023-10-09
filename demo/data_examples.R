
source("my_ets.R")

library(doParallel)
library(tidyverse)

error_measures <- function(a,f){
  return(c(
    smape = 2 * mean(abs(a - f) / (abs(a) + abs(f))),
    mae = mean(abs(a - f)),
    rmse = sqrt(mean((a - f) ^ 2))
  ))
}

# data preparation ----

# M3
#install.packages("Mcomp")

m3_type <- c("DEMOGRAPHIC", "FINANCE", "INDUSTRY", "MACRO", "MICRO")
datasets <- sapply(m3_type, function(t) {
  d <- subset(Mcomp::M3, "MONTHLY", t)
  all_dates <-
    lapply(d, function(ts)
      seq(as.Date(paste(
        c(start(ts$x), "01"), collapse = "-"
      )), as.Date(paste(
        c(end(ts$xx), "01"), collapse = "-"
      )), by = "1 month"))
  all_dates <- sort(unique(do.call(c, all_dates)))


  data <- lapply(d, function(ts) {
    dates <-
      seq(as.Date(paste(c(start(
        ts$x
      ), "01"), collapse = "-")), as.Date(paste(c(end(
        ts$xx
      ), "01"), collapse = "-")), by = "1 month")

    cbind(which(all_dates %in% dates), c(ts$x, ts$xx))
  })

  test_idx <- unlist(lapply(d, function(ts) {
    which(all_dates == as.Date(paste(c(
      start(ts$xx), "01"
    ), collapse = "-")))
  }))

  tibble(data,test_idx,min_ts=round(0.8*length(all_dates)))

}, simplify=FALSE)

# Tourism Data
# install.packages("Tcomp")
tourism <- subset(Tcomp::tourism,"monthly")
all_dates <-
  lapply(tourism, function(ts)
    seq(as.Date(paste(
      c(start(ts$x), "01"), collapse = "-"
    )), as.Date(paste(
      c(end(ts$xx), "01"), collapse = "-"
    )), by = "1 month"))
all_dates <- sort(unique(do.call(c, all_dates)))

data <- lapply(tourism, function(ts) {
  dates <-
    seq(as.Date(paste(c(
      start(ts$x), "01"
    ), collapse = "-")), as.Date(paste(c(
      end(ts$xx), "01"
    ), collapse = "-")), by = "1 month")

  cbind(which(all_dates %in% dates), c(ts$x, ts$xx))
})

test_idx <- unlist(lapply(tourism, function(ts) {
  which(all_dates == as.Date(paste(c(
    start(ts$xx), "01"
  ), collapse = "-")))
}))

datasets <- c(datasets, list(TOURISM = tibble(data, test_idx, min_ts=0.8*length(all_dates))))

# CIF 2016
cif2016 <- readRDS("data/cif2016.rds")
cif2016_processed <- cif2016[[1]] %>%
  group_by(series_name) %>%
  mutate(length = n()) %>%
  ungroup() %>%
  mutate(max_length=max(length)) %>%
  group_by(series_name) %>%
  reframe(data=list(cbind((max_length[1]-length(series_value)+1):max_length[1], series_value)), # the ts end at the same time?
            test_idx=max_length[1]-as.numeric(horizon)+1,
          min_ts=round(0.8*max_length[1])) %>%
  distinct()

datasets <- c(datasets, list(CIF2016 = cif2016_processed))

# Hospital Data
# install.packages("expsmooth")

hospital <- tibble(data=apply(expsmooth::hospital,2,function(col)cbind(1:nrow(expsmooth::hospital),col,deparse.level = 0),simplify=FALSE),
                   test_idx = 73,
                   min_ts = round(0.8*84))

datasets <- c(datasets, list(HOSPITAL = hospital))


# experiment ----
result <- list()
i<-1
model <- list(model = my_ets,
              forecast = forecast.my_ets,
              refit = refit.my_ets,
              fitted = fitted.my_ets,
              residuals = residuals.my_ets)
for (ds in datasets) {
  res_cv  <-
    ts_avg.cv(
      ds$data,
      model,
      xtest_idx = ds$test_idx,
      min_ts = unique(ds$min_ts),
      type = c("simple_avg","simple_avg_n", "dist_avg_n", "dist_avg","global_avg","perf_avg","perf_avg_r"),
      c(1, 5, 10, 20),
      test_multistep = TRUE,
      parallel = 6,
      benchmark = "model"
    )

  res_best_onestep <- ts_avg(
    ds$data,
    model,
    k = res_cv$final_model$k,
    type = res_cv$final_model$type,
    xtest_idx = ds$test_idx,
    benchmark = "model",
    test_multistep = FALSE
  )

  eval_multistep <- ds %>%
    ungroup() %>%
    mutate(tsavg = lapply(res_cv$final_model$forecasts,function(f)f[1,-dim(f)[2],2])) %>%
    rowwise() %>%
    mutate(ets = list({
      m_ets <- forecast::ets(data[data[,1]<test_idx,2])
      as.numeric(forecast::forecast(m_ets,h=max(data[,1])-test_idx+1)$mean)
    }),
    arima = list({
      m_arima <- forecast::auto.arima(data[data[,1]<test_idx,2])
      as.numeric(forecast::forecast(m_arima,h=max(data[,1])-test_idx+1)$mean)
    })) %>%
    pivot_longer(c(tsavg, ets, arima), names_to="model") %>%
    rowwise() %>%
    mutate(err=list(error_measures(data[data[,1]>=test_idx,2],as.numeric(value)))) %>%
    unnest_wider(err)

  eval_onestep <- ds %>%
    ungroup() %>%
    mutate(tsavg = lapply(res_best_onestep$forecasts,function(f)f[-dim(f)[1],1,2])) %>%
    rowwise() %>%
    mutate(ets = list({
      sapply(test_idx:max(data[,1]),function(ind){
        m_ets <- forecast::ets(data[data[,1]<ind,2])
        as.numeric(forecast::forecast(m_ets,h=1)$mean)
      })
    }),
    arima = list({
      sapply(test_idx:max(data[, 1]), function(ind) {
        m_arima <- forecast::auto.arima(data[data[, 1] < ind, 2])
        as.numeric(forecast::forecast(m_arima, h = 1)$mean)
      })
    })) %>%
    pivot_longer(c(tsavg, ets, arima), names_to="model") %>%
    rowwise() %>%
    mutate(err=list(error_measures(data[data[,1]>=test_idx,2],as.numeric(value)))) %>%
    unnest_wider(err)


  result <- c(result, list(res_multistep = res_cv$final_model,
                           res_onestep = res_best_onestep,
                           eval_multistep = eval_multistep,
                           eval_onestep = eval_onestep))
  message(sprintf("subset %d/%d finished", i,length(datasets)))
  i<-i+1
}

names(result) <- names(datasets)

saveRDS(result,file="./result_datasets.rds")
