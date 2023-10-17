simple_avg_fn <- function(x, models, neighbors){
  refit_models <- lapply(models$models[-1], models$refit_fn, ynew=x[[1]][,-1])
  forecasts <- c(list(models$forecast_fn(models$models[[1]], h=models$h)),
                 lapply(refit_models,models$forecast_fn,h=models$h))

  return(Reduce('+', forecasts)/length(forecasts))
}

simple_avg_n_fn <- function(x, models, neighbors){
  refit_models <- lapply(models$models[-1], models$refit_fn, ynew=x[[1]][,-1])
  forecasts <- lapply(refit_models,models$forecast_fn,h=models$h)

  return(Reduce('+', forecasts)/length(forecasts))
}

dist_avg_n_fn <- function(x, models, neighbors){
  refit_models <- lapply(models$models[-1], models$refit_fn, ynew=x[[1]][,-1])
  forecasts <- lapply(refit_models,models$forecast_fn,h=models$h)
  weights <- 1/neighbors$dists
  weights <- weights/sum(weights)

  return(Reduce('+', Map('*', forecasts, weights)))
}

dist_avg_fn <- function(x, models, neighbors){
  refit_models <- lapply(models$models[-1], models$refit_fn, ynew=x[[1]][,-1])
  forecasts <- c(list(models$forecast_fn(models$models[[1]],h=models$h)),
                 lapply(refit_models,models$forecast_fn,h=models$h))

  init_cent <- neighbors$data[[which.max(lapply(neighbors$data,nrow))]]
  avg <- adba(neighbors$data, init_cent)

  weights <- 1/avg$dist_to_avg
  weights <- weights/sum(weights)

  return(Reduce('+', Map('*', forecasts, weights)))
}

global_avg_fn <- function(x, models, neighbors){
  init_cent <- neighbors$data[[which.max(lapply(neighbors$data,nrow))]]
  avg <- adba(neighbors$data, init_cent)
  avg_model <- models$model_fn(avg$avg[,-1])

  return(models$forecast_fn(models$refit_fn(avg_model,x[[1]][,-1]), h=models$h))
}

perf_avg_fn <- function(x, models, neighbors){
  refit_models <- lapply(models$models[-1], models$refit_fn, ynew=x[[1]][,-1])
  forecasts <- c(list(models$forecast_fn(models$models[[1]],h=models$h)),
                 lapply(refit_models,models$forecast_fn,h=models$h))

  refit_models_2 <- sapply(seq_along(x), function(i){
    models$refit_fn(models$models[[i]], ynew=x[[i]][,-1])
  },simplify=FALSE)

  weights <- sapply(seq_along(x), function(i) {
    yhat <- cbind(x[[i]][,1],models$fitted_fn(refit_models_2[[i]],h=models$h))
    1/rmsse(x[[i]], yhat)
  })

  weights <- weights/sum(weights)

  return(Reduce('+', Map('*', forecasts, weights)))

}

perf_avg_r_fn <- function(x, models, neighbors){
  refit_models <- lapply(models$models, models$refit_fn, ynew=x[[1]][,-1])
  forecasts <- lapply(refit_models,models$forecast_fn,h=models$h)

  weights <- sapply(seq_along(x), function(i) {
    fitted_vals <- models$fitted_fn(refit_models[[i]],h=models$h)
    if (all(is.na(fitted_vals))) fitted_vals <- matrix(NA, nrow=nrow(x[[1]]), ncol=ncol(x[[1]])-1)
    yhat <- cbind(x[[1]][,1], fitted_vals)
    1/rmsse(x[[1]], yhat)
  })

  weights <- weights/sum(weights)

  return(Reduce('+', Map('*', forecasts, weights)))

}
