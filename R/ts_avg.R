#' ts_avg
#'
#' Perform TS Avg Forecasting for a particular type
#' @param x list of matrices. each matrix corresponds to one time series whereby the first column always corresponds to its index.
#' @param model_list a generic model function. the resulting model object must have forecast, refit, fitted, and residuals functions (as demonstrated in my_ets).
#' @param k vector of number of neighbors used for each time series in x.
#' @param h forecast horizon. if h>1, then multiple forecasts are computed. default=1.
#' @param type averaging type. one of  "simple_avg", "simple_avg_n", "dist_avg_n", "dist_avg", "global_avg", "perf_avg", "perf_avg_r".
#' @param xtest_idx index to split each time series in train and test set. default = NULL.
#' @param test_multistep if TRUE, then multistep-ahead forecasts are performed on test set.
#' @param benchmark benchmark forecast type for computing the errors. either random walk forecast "rw" or model forecast "model".
#' @param detailed_result if TRUE, a more detailed result object is returned.
#' @param execution_time if TRUE, also return execution time.
#' @param ... additional parameters for model function
#' @return object of class ts_avg containing
#' \itemize{
#'  \item data
#'  \item type
#'  \item forecasts
#'  \item h
#'  \item test_multistep
#'  \item xtest_idx
#'  \item test_errors
#'  \item cum_test_errors
#'  \item neighbors (if detailed_result is TRUE)
#'  \item exec_time (if execution_time is TRUE)
#' }
#' @example demo/demo.R
#' @seealso [plot.ts_avg]
#' @export

ts_avg <- function(x,
                   model_list,
                   k,
                   h = 1,
                   type=c("simple_avg", "simple_avg_n", "dist_avg_n", "dist_avg","global_avg","perf_avg","perf_avg_r"),
                   xtest_idx=NULL,
                   test_multistep=FALSE,
                   benchmark=c("rw","model"),
                   detailed_result=FALSE,
                   execution_time=FALSE,
                   ...){
  stopifnot(all(lapply(x,ncol)>1))

  type <- match.arg(type)
  benchmark <- match.arg(benchmark)

  avg_fn <- get(paste(type,"fn",sep="_"), envir = getNamespace("TSAvg"))

  if (execution_time){
    start_time <- Sys.time()
  }

  if (!is.null(xtest_idx)){

    if (length(xtest_idx) == 1){
      xtest_idx <- rep(xtest_idx, length(x))
    }

    xtest <- sapply(seq_along(x), function(i) x[[i]][x[[i]][,1] >= xtest_idx[i], ,drop=FALSE], simplify=FALSE)
    xtrain <- sapply(seq_along(x), function(i) x[[i]][x[[i]][,1] < xtest_idx[i], ,drop=FALSE], simplify=FALSE)
    test_idx <- lapply(xtest, function(xti){
      if (nrow(xti)==0) return(NULL)
      else return(xti[,1])
    })
  }
  else {
    xtrain <- x
    test_idx <- replicate(length(x), NULL, simplify = FALSE)
  }

  models <- lapply(xtrain, function(xi) model_list$model(xi[,-1], ...))

  if (length(k)==1) k <- rep(k,length(x))

  neighbors <- dtw_neighbors(xtrain, k)

  result_fc <- list()
  result_fitted <- list()

  for (i in seq_along(x)) {
    if (test_multistep == FALSE) {
      result_i <-
        array(dim = c(length(test_idx[[i]]) + 1, h, ncol(xtrain[[i]])))

      fitted_i <- array(dim = c(nrow(xtrain[[i]]), ncol(xtrain[[i]])))
      ndata <- lapply(neighbors, '[[', i)
      if (k[i] > 0 & any(!is.infinite(ndata$dists))) {
        ndata$data <- neighbors$data[c(i, ndata$neighbors)]
        tmp_res <-
          avg_fn(c(xtrain[i], xtrain[neighbors$neighbors[[i]]]),
                 list(
                   model_fn = function(x)
                     model_list$model(x, ...),
                   refit_fn = model_list$refit,
                   forecast_fn = model_list$forecast,
                   fitted_fn = model_list$fitted,
                   models = c(models[i], models[neighbors$neighbors[[i]]]),
                   h = h
                 ),
                 ndata)
      }
      else {
        tmp_res <- list(forecasts=model_list$forecast(models[[i]], h = h),
                        fitted_values=model_list$fitted(models[[i]], h = h))
      }

      fitted_i <- tmp_res$fitted_values

      result_i[1, , ] <-
        cbind(max(xtrain[[i]][, 1]) + (1:h), tmp_res$forecasts)

      if (!is.null(xtest_idx)) {
        for (j in seq_along(test_idx[[i]])) {
          tj <- test_idx[[i]][j]
          x2 <-
            list(rbind(xtrain[[i]], xtest[[i]][xtest[[i]][, 1] <= tj, ]))
          x2_neighbors <- sapply(neighbors$neighbors[[i]],
                                 function(n)
                                   rbind(xtrain[[n]], xtest[[n]][xtest[[n]][, 1] <= tj, ]), simplify = FALSE)
          if (k[i] > 0) {
            ndata <- lapply(neighbors, '[[', i)
            ndata$data <- Map(function(y, c) {
              y[, -1] <- y[, -1] - c
              return(y)
            },
            c(x2, x2_neighbors),
            neighbors$centers[c(i, ndata$neighbors)])
            tmp_res <- avg_fn(
              c(x2, x2_neighbors),
              list(
                model_fn = function(x)
                  model_list$model(x, ...),
                refit_fn = model_list$refit,
                forecast_fn = model_list$forecast,
                fitted_fn = model_list$fitted,
                models = c(models[i], models[neighbors$neighbors[[i]]]),
                h = h
              ),
              ndata
            )
          }
          else {
            tmp_res <- list(forecasts=model_list$forecast(model_list$refit(models[[i]], x2[[1]][, -1]), h = h))
          }
          result_i[j + 1, , ] <-
            cbind(max(x2[[1]][, 1]) + (1:h), tmp_res$forecasts)
        }
      }
    }
    else {
      h_i <- length(test_idx[[i]]) + 1
      result_i <- array(dim = c(1, h_i, ncol(xtrain[[i]])))

      ndata <- lapply(neighbors, '[[', i)
      if (k[i] > 0 & any(!is.infinite(ndata$dists))) {
        ndata$data <- neighbors$data[c(i, ndata$neighbors)]
        tmp_res <-
          avg_fn(c(xtrain[i], xtrain[neighbors$neighbors[[i]]]),
                 list(
                   model_fn = function(x)
                     model_list$model(x, ...),
                   refit_fn = model_list$refit,
                   forecast_fn = model_list$forecast,
                   fitted_fn = model_list$fitted,
                   models = c(models[i], models[neighbors$neighbors[[i]]]),
                   h = h_i
                 ),
                 ndata)

      }
      else {
        tmp_res <- list(forecasts=model_list$forecast(models[[i]], h = h_i))
      }

      result_i[1, , ] <- cbind(max(xtrain[[i]][, 1]) + (1:h_i), tmp_res$forecasts)
    }
    result_fc <- c(result_fc, list(result_i))
    result_fitted <- c(result_fitted, list(fitted_i))
  }

  if (is.null(xtest_idx)){
    #result <- lapply(result, "[",1,TRUE,TRUE)
    result_fc <- lapply(result_fc, function(r){
      matrix(r[1, , ],ncol=dim(r)[3])
    })
    cum_test_errors <- test_errors <- NULL
  }
  else {
    agg_dim <- ifelse(test_multistep, 1, 2)
    cum_test_errors <- sapply(seq_along(result_fc), function(i) {
      apply(result_fc[[i]], agg_dim, function(r) {
        ind <- intersect(r[,1], xtest[[i]][,1])
        if (length(ind) == 0) return(NA)
        scaling <- if (benchmark == "rw")
          rmse(xtrain[[i]], rw(xtrain[[i]]))
        else
          rmse(cbind(0, model_list$residuals(models[[i]])), 0)
        rmsse(xtest[[i]][xtest[[i]][,1] %in% ind, , drop=FALSE], r[r[,1] %in% ind, , drop=FALSE],
              scaling = scaling,
              cumulative = TRUE)
      })
    }, simplify=FALSE)

    test_errors <- sapply(seq_along(result_fc), function(i) {
      apply(result_fc[[i]], agg_dim, function(r) {
        ind <- intersect(r[,1], xtest[[i]][,1])
        if (length(ind) == 0) return(NA)
        scaling <- if (benchmark == "rw")
          rmse(xtrain[[i]], rw(xtrain[[i]]))
        else
          rmse(cbind(0, model_list$residuals(models[[i]])), 0)
        rmsse(xtest[[i]][xtest[[i]][,1] %in% ind, , drop=FALSE], r[r[,1] %in% ind, , drop=FALSE],
              scaling = scaling,
              cumulative = FALSE)
      })
    }, simplify=FALSE)
  }

  if (execution_time){
    time_needed <- as.numeric(Sys.time() - start_time,"secs")
  }

  res <- list(data = x,
              type = type,
              k=k,
              forecasts = result_fc,
              fitted = result_fitted,
              h=h,
              test_multistep=test_multistep,
              xtest_idx=xtest_idx,
              test_errors = test_errors,
              cum_test_errors = cum_test_errors)

  if (detailed_result){
    res <- c(res, list(neighbors = neighbors))
  }

  if (execution_time){
    res <- c(res, list(exec_time = time_needed))
  }

  attr(res,"class") <- "ts_avg"

  return(res)
}

#' plot.ts_avg
#'
#' Plot ts_avg result
#' @param object result of ts_avg function
#' @param plot_type response or error
#' @param h forecast horizon to use. default=1.
#' @param ids vector of ids to filter appropriate plots.
#' @param ref_plot ts_avg plots can be concatenated such that different types and ks can be shwon together.
#' @return ggplot2 object
#' @example demo/demo.R
#' @import ggplot2
#' @import reshape2
#' @export

plot.ts_avg <- function(object, plot_type=c("response","error","error_quants"), h=1, ids=NULL, ref_plot=NULL){

  plot_type <- match.arg(plot_type)
  if (plot_type == "response") {
    object$data <- lapply(object$data, function(d){
      colnames(d) <- NULL
      return(d)
    })
    x_df <- cbind(avg_type=NA,
                  forecast = FALSE,
                  as.data.frame(do.call(
                    rbind, sapply(seq_along(object$data), function(i)
                      cbind(id = i, fc_idx = NA, object$data[[i]]),
                      simplify=FALSE)
                  )))

    if (object$h == 1) {
      warning(sprintf("parameter h=%d is not used", h))
      h_plot <- ifelse(object$test_multistep, "multistep", 1)
    }
    else {
      h_plot <- h
    }


    fc_df <- if (object$test_multistep == FALSE){
      cbind(avg_type=object$type,
                   forecast = TRUE,
                   as.data.frame(do.call(
                     rbind, sapply(seq_along(object$data), function(i) {
                       fc <- if (is.null(object$xtest_idx)) {
                         matrix(object$forecasts[[i]], ncol = ncol(object$forecasts[[i]]))
                       }
                       else {
                         matrix(object$forecasts[[i]][, h_plot, ], ncol=dim(object$forecasts[[i]])[3])
                       }
                       cbind(id = i, fc_idx = min(fc[, 1]) - h_plot, fc)
                     }, simplify = FALSE)
                   )))
    }
    else {
      cbind(avg_type=object$type,
            forecast = TRUE,
            as.data.frame(do.call(
              rbind, sapply(seq_along(object$data), function(i) {
                fc <- if (is.null(object$xtest_idx)) {
                  matrix(object$forecasts[[i]], ncol = ncol(object$forecasts[[i]]))
                }
                else {
                  matrix(object$forecasts[[i]][1, , ], ncol=dim(object$forecasts[[i]])[3])
                }
                cbind(id = i, fc_idx = min(fc[, 1]) - 1, fc)
              }, simplify = FALSE)
            )))
    }

    plot_df <- rbind(x_df, fc_df)

    if (!is.null(ids)){
      plot_df <- plot_df[plot_df$id %in% ids, ]
    }

    colnames(plot_df)[1:5] <- c("avg_type", "forecast", "id", "fc_idx", "index")
    colnames(plot_df)[-(1:5)] <- paste0("V", 1:(ncol(plot_df) - 5))

    #plot_res(plot_df, h_plot)
    get("plot_res",envir = getNamespace("TSAvg"))(plot_df, h_plot)
  }
  else if (plot_type == "error"){
    error_df <- as.data.frame(cbind(1:length(object$test_errors), do.call(rbind,object$test_errors)[,h]))
    #error_df <- error_df[order(error_df[,2]),]
    colnames(error_df) <- c("id","rmsse")
    error_df$id <- factor(error_df$id, levels=error_df[order(error_df[,2]),1])

    if (!is.null(ids)){
      error_df <- error_df[error_df$id %in% ids, ]
    }

    ggplot2::ggplot(error_df, ggplot2::aes(x=id,y=rmsse, group=1))+
      ggplot2::geom_step()+
      ggplot2::geom_hline(yintercept=1,linetype="dashed")
  }
  else if (plot_type == "error_quants"){
    error_df <- as.data.frame(cbind(1:length(object$test_errors), do.call(rbind,object$test_errors)[,h]))
    #error_df <- error_df[order(error_df[,2]),]
    colnames(error_df) <- c("id","rmsse")
    error_df$id <- factor(error_df$id, levels=error_df[order(error_df[,2]),1])
    error_df$q <- rank(error_df[,2])/nrow(error_df)
    error_df$type <- object$type
    error_df$mean_k <- mean(object$k)
    # if (!is.null(ids)){
    #   error_df <- error_df[error_df$id %in% ids, ]
    # }

    if (is.null(ref_plot)){
      p <-
        ggplot2::ggplot(error_df,
                        ggplot2::aes(
                          x = q,
                          y = rmsse,
                          #group = 1,
                          color = interaction(type, round(mean_k, 1),sep=", ")
                        )) +
        ggplot2::geom_step() +
        ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
        ggplot2::geom_vline(ggplot2::aes(xintercept = max(q[rmsse < 1]))) +
        ggplot2::scale_x_continuous(name = "Percentiles", breaks = (1:10) / 10) +
        ggplot2::scale_color_discrete("Type, Avg. k")

    }
    else {
      error_df <- rbind(p$data, error_df)
      p <-
        ggplot2::ggplot(error_df,
                        ggplot2::aes(
                          x = q,
                          y = rmsse,
                          #group = 1,
                          color = interaction(type, round(mean_k, 1), sep=", ")
                        )) +
        ggplot2::geom_step() +
        ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
        ggplot2::geom_vline(data = error_df, ggplot2::aes(xintercept = max(q[rmsse < 1]))) +
        ggplot2::scale_x_continuous(name = "Percentiles", breaks = (1:10) / 10)+
        ggplot2::scale_color_discrete("Type, Avg. k")
    }
    return(p)
  }
}
