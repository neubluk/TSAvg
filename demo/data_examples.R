library(doParallel)
library(tidyverse)
library(TSAvg)

error_measures <- function(a,f,rw_train_error=NULL){

  if (is.null(rw_train_error)){
    scalings <- sapply(seq_along(a), function(i) sqrt(mean( (a[1:i]-lag(a[1:i]))^2, na.rm=TRUE)))
    rmsse <- sqrt(mean( ((a-f)/scalings)^2, na.rm=TRUE))
  }
  else {
    rmsse <- sqrt(mean( ((a-f)/rw_train_error)^2, na.rm=TRUE))
  }

  return(c(
    smape = 2 * mean(abs(a - f) / (abs(a) + abs(f)), na.rm=TRUE),
    mae = mean(abs(a - f), na.rm=TRUE),
    rmse = sqrt(mean((a - f) ^ 2, na.rm=TRUE)),
    rmsse = rmsse
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

  tibble(id = 1:length(data), data,test_idx,min_ts=round(0.8*min(test_idx)), max_input = round(1.25*max(max(unlist(lapply(d,function(ts)ts$h))),12)))

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

datasets <- c(datasets, list(TOURISM = tibble(id=1:length(data),data, test_idx, min_ts=round(0.8*min(test_idx)),max_input = round(1.25*max(max(unlist(lapply(tourism,function(ts)ts$h))),12)))))

# CIF 2016
data("cif2016")
cif2016_processed <- cif2016[[1]] %>%
  group_by(series_name) %>%
  mutate(length = n()) %>%
  ungroup() %>%
  mutate(max_length=max(length)) %>%
  group_by(series_name) %>%
  reframe(data=list(cbind((max_length[1]-length(series_value)+1):max_length[1], series_value)), # the ts end at the same time?
            test_idx=max_length[1]-as.numeric(horizon)+1,
          min_ts=round(0.8*max_length[1]),
          max_input = round(1.25*max(as.numeric(horizon),12))) %>%
  distinct() %>%
  rename(id=series_name)

datasets <- c(datasets, list(CIF2016 = cif2016_processed))

# Hospital Data
# install.packages("expsmooth")

hospital <- tibble(data=apply(expsmooth::hospital,2,function(col)cbind(1:nrow(expsmooth::hospital),col,deparse.level = 0),simplify=FALSE),
                   test_idx = 73,
                   min_ts = round(0.8*73),
                   max_input = round(1.25*12)) %>%
  mutate(id=1:length(data), .before=data)

datasets <- c(datasets, list(HOSPITAL = hospital))

# Food Demand Data

data("food_demand")
food_demand_data <- food_demand %>%
  group_by(fridge_id) %>%
  summarise(data=list(cbind(index,sold)),
            min_ts=27,
            test_idx=81,
            max_input = 5) %>%
  rename(id=fridge_id)

datasets <- c(datasets, list(FOOD_DEMAND = food_demand_data))

# experiment ----
source("demo/my_ets.R")
model <- list(model = my_ets,
              forecast = forecast.my_ets,
              refit = refit.my_ets,
              fitted = fitted.my_ets,
              residuals = residuals.my_ets)
types <- c(
  "simple_avg",
  "simple_avg_n",
  "dist_avg_n",
  "dist_avg",
  "global_avg",
  "perf_avg",
  "perf_avg_r"
)

nr_cores <- 4
result <- list()
i <- 1
for (ds in datasets) {
  lambda <- NULL
  #lambda <- if (min(unlist(lapply(ds$data,function(d)min(d[,2]))))>0) 0 # need to run this again with log trafos if applicable!!!
  job::job({
    res_cv  <-
      ts_avg.cv(
        ds$data,
        model,
        xtest_idx = ds$test_idx,
        min_ts = unique(ds$min_ts),
        type = types,
        c(1, 5, 10, 20),
        test_multistep = FALSE,
        parallel = nr_cores,
        benchmark = "model",
        lambda = lambda
      )

    best_result_onestep <- lapply(types, function(t) {
      ts_avg(
        ds$data,
        model,
        res_cv$best_k[, t],
        h = 1,
        type = t,
        xtest_idx = ds$test_idx,
        benchmark = "model",
        test_multistep = TRUE
      )
    })

    names(best_result_onestep) <- types

    best_type_onestep <- which.min(lapply(best_result_onestep,function(br)mean(unlist(br$test_errors))))
    final_model_onestep <- best_result_onestep[[best_type_onestep]]

    eval_multistep <- ds %>%
      ungroup() %>%
      mutate(tsavg = lapply(res_cv$final_model$forecasts, function(f)
        f[1, -dim(f)[2], 2])) %>%
      rowwise() %>%
      mutate(ets = list({
        m_ets <- forecast::ets(data[data[, 1] < test_idx, 2])
        as.numeric(forecast::forecast(m_ets, h = max(data[, 1]) - test_idx +
                                        1)$mean)
      }),
      arima = list({
        m_arima <- forecast::auto.arima(data[data[, 1] < test_idx, 2])
        as.numeric(forecast::forecast(m_arima, h = max(data[, 1]) - test_idx +
                                        1)$mean)
      })) %>%
      pivot_longer(c(tsavg, ets, arima), names_to = "model") %>%
      rowwise() %>%
      mutate(err = list(error_measures(data[data[, 1] >= test_idx, 2], as.numeric(value)))) %>%
      unnest_wider(err)

    eval_onestep <- ds %>%
      ungroup() %>%
      mutate(tsavg = lapply(final_model_onestep$forecasts, function(f) # change back to final_model_onestep
        f[-dim(f)[1], 1, 2])) %>%
      rowwise() %>%
      mutate(ets = list({
        sapply(test_idx:max(data[, 1]), function(ind) {
          m_ets <- forecast::ets(data[data[, 1] < ind, 2], lambda=lambda)
          as.numeric(forecast::forecast(m_ets, h = 1)$mean)
        })
      }),
      arima = list({
        sapply(test_idx:max(data[, 1]), function(ind) {
          m_arima <- forecast::auto.arima(data[data[, 1] < ind, 2], lambda = lambda)
          as.numeric(forecast::forecast(m_arima, h = 1)$mean)
        })
      })) %>%
      pivot_longer(c(tsavg, ets, arima), names_to = "model") %>%
      rowwise() %>%
      mutate(err = list(error_measures(data[data[, 1] >= test_idx, 2], as.numeric(value)))) %>%
      unnest_wider(err)

    saveRDS(
      list(
        #res_onestep = res_cv,
        res_multistep = res_cv,
        res_onestep = final_model_onestep,
        eval_multistep = eval_multistep,
        eval_onestep = eval_onestep
      ),
      file = sprintf("./test_%s.rds", names(datasets)[i])
    )

    job::export("none")
  },
  title=names(datasets)[i],
  import="all")
  message(sprintf("subset %d/%d finished", i,length(datasets)))
  i<-i+1
}



# ids <- c(3,14,25,27)

res_cv  <-
  ts_avg.cv(
    food_demand_data$data,
    model,
    xtest_idx = food_demand_data$test_idx,
    min_ts = unique(food_demand_data$min_ts),
    type = c(
      "simple_avg",
      "simple_avg_n",
      "dist_avg_n",
      "dist_avg",
      "global_avg",
      "perf_avg",
      "perf_avg_r"
    ),
    c(1, 5, 10, 20),
    test_multistep = FALSE,
    parallel = FALSE,
    benchmark = "model",
    lambda=0
  )


