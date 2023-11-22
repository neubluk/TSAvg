simple_avg_fn <- function(x, models, neighbors){
  refit_models <- lapply(models$models, models$refit_fn, ynew=x[[1]][,-1])
  forecasts <- lapply(refit_models, models$forecast_fn, h = models$h)
  fitted_values <- lapply(refit_models, models$fitted_fn, h = models$h)

  return(list(forecasts=Reduce('+', forecasts)/length(forecasts),
              fitted_values=Reduce('+', fitted_values)/length(fitted_values)))
}

simple_avg_n_fn <- function(x, models, neighbors){
  refit_models <- lapply(models$models[-1], models$refit_fn, ynew=x[[1]][,-1])
  forecasts <- lapply(refit_models,models$forecast_fn,h=models$h)
  fitted_values <- lapply(refit_models, models$fitted_fn, h = models$h)

  return(list(forecasts=Reduce('+', forecasts)/length(forecasts),
              fitted_values=Reduce('+', fitted_values)/length(fitted_values)))
}

dist_avg_n_fn <- function(x, models, neighbors){
  refit_models <- lapply(models$models[-1], models$refit_fn, ynew=x[[1]][,-1])
  forecasts <- lapply(refit_models,models$forecast_fn,h=models$h)
  fitted_values <- lapply(refit_models,models$fitted_fn,h=models$h)
  weights <- 1/neighbors$dists
  weights <- weights/sum(weights)

  return(list(forecasts=Reduce('+', Map('*', forecasts, weights)),
              fitted_values=Reduce('+', Map('*', fitted_values, weights))))
}

dist_avg_fn <- function(x, models, neighbors){
  refit_models <- lapply(models$models, models$refit_fn, ynew=x[[1]][,-1])
  forecasts <- lapply(refit_models,models$forecast_fn,h=models$h)

  fitted_values <- lapply(refit_models,models$fitted_fn,h=models$h)

  init_cent <- neighbors$data[[which.max(lapply(neighbors$data,nrow))]]
  avg <- adba(neighbors$data, init_cent)

  weights <- 1/avg$dist_to_avg
  weights <- weights/sum(weights)

  return(list(forecasts=Reduce('+', Map('*', forecasts, weights)),
              fitted_values=Reduce('+', Map('*', fitted_values, weights))))
}

global_avg_fn <- function(x, models, neighbors){
  init_cent <- neighbors$data[[which.max(lapply(neighbors$data,nrow))]]
  avg <- adba(neighbors$data, init_cent)
  avg_model <- models$model_fn(avg$avg[,-1])

  avg_model_refit <- models$refit_fn(avg_model,x[[1]][,-1])

  return(list(forecasts=models$forecast_fn(avg_model_refit, h=models$h),
              fitted_values=models$fitted_fn(avg_model_refit, h=models$h)))
}

perf_avg_fn <- function(x, models, neighbors){
  refit_models <- lapply(models$models, models$refit_fn, ynew=x[[1]][,-1])
  forecasts <- lapply(refit_models,models$forecast_fn,h=models$h)

  fitted_values <- lapply(refit_models,models$fitted_fn,h=models$h)

  refit_models_2 <- sapply(seq_along(x), function(i){
    models$refit_fn(models$models[[i]], ynew=x[[i]][,-1])
  },simplify=FALSE)

  weights <- sapply(seq_along(x), function(i) {
    yhat <- cbind(x[[i]][,1],models$fitted_fn(refit_models_2[[i]],h=models$h))
    1/rmsse(x[[i]], yhat)
  })

  weights <- weights/sum(weights)

  return(list(forecasts=Reduce('+', Map('*', forecasts, weights)),
              fitted_values=Reduce('+', Map('*', fitted_values, weights))))

}

perf_avg_r_fn <- function(x, models, neighbors){
  refit_models <- lapply(models$models, models$refit_fn, ynew=x[[1]][,-1])
  forecasts <- lapply(refit_models,models$forecast_fn,h=models$h)
  fitted_values <- lapply(refit_models,models$fitted_fn,h=models$h)

  weights <- sapply(seq_along(x), function(i) {
    fitted_vals <- fitted_values[[i]]
    if (all(is.na(fitted_vals))) fitted_vals <- matrix(NA, nrow=nrow(x[[1]]), ncol=ncol(x[[1]])-1)
    yhat <- cbind(x[[1]][,1], fitted_vals)
    1/rmsse(x[[1]], yhat)
  })

  weights <- weights/sum(weights)

  return(list(forecasts=Reduce('+', Map('*', forecasts, weights)),
              fitted_values=Reduce('+', Map('*', fitted_values, weights))))

}
