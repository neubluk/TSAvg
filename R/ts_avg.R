#' ts_avg
#'
#' Perform TS Avg Forecasting for a particular type
#' @param x list of matrices. each matrix corresponds to one time series whereby the first column always corresponds to its index.
#' @param model a generic model function. the resulting model object must have forecast, refit, fitted, and residuals functions (as demonstrated in my_ets).
#' @param k vector of number of neighbors used for each time series in x.
#' @param h forecast horizon. if h>1, then multiple forecasts are computed. default=1.
#' @param type averaging type. one of  "simple_avg", "simple_avg_n", "dist_avg_n", "dist_avg", "global_avg", "perf_avg", "perf_avg_r".
#' @param xtest_idx index to split the every time series in train and test set. default = NULL.
#' @param benchmark benchmark forecast type for computing the errors. either random walk forecast "rw" or model forecast "model".
#' @param ... additional parameters for model function
#' @return object of class ts_avg
#' @example demo/demo.R
#' @seealso [plot.ts_avg]
#' @export

ts_avg <- function(x,
                   model,
                   k,
                   h = 1,
                   type=c("simple_avg", "simple_avg_n", "dist_avg_n", "dist_avg","global_avg","perf_avg","perf_avg_r"),
                   xtest_idx=NULL,
                   benchmark=c("rw","model"),
                   ...){
  stopifnot(all(lapply(x,ncol)>1))

  type <- match.arg(type)
  benchmark <- match.arg(benchmark)

  avg_fn <- get(paste(type,"fn",sep="_"), envir = getNamespace("TSAvg"))

  if (!is.null(xtest_idx)){
    xtest <- lapply(x, function(xi)xi[xi[,1] >= xtest_idx, ,drop=FALSE])
    xtrain <- lapply(x, function(xi)xi[xi[,1] < xtest_idx, ,drop=FALSE])
    test_idx <- lapply(xtest, function(xti){
      if (nrow(xti)==0) return(NULL)
      else return(xti[,1])
    })
  }
  else {
    xtrain <- x
    test_idx <- replicate(length(x), NULL, simplify = FALSE)
  }

  models <- lapply(xtrain, function(xi) model(xi[,-1], ...))

  if (length(k)==1) k <- rep(k,length(x))

  neighbors <- dtw_neighbors(xtrain, k)

  result <- list()

  for (i in seq_along(x)){
    result_i <- array(dim=c(length(test_idx[[i]])+1,h,ncol(xtrain[[i]])))
    ndata <- lapply(neighbors, '[[', i)
    if (k[i] > 0 & any(!is.infinite(ndata$dists))) {
      ndata$data <- neighbors$data[c(i,ndata$neighbors)]
      tmp_res <- avg_fn(
        c(xtrain[i], xtrain[neighbors$neighbors[[i]]]),
        list(
          model_fn = function(x)
            model(x, ...),
          models = c(models[i], models[neighbors$neighbors[[i]]]),
          h = h
        ),
        ndata
      )

    }
    else {
      tmp_res <- forecast(models[[i]], h = h)
    }

    result_i[1, , ] <- cbind(max(xtrain[[i]][,1])+(1:h), tmp_res)

    if (!is.null(xtest_idx)) {
      for (j in seq_along(test_idx[[i]])){
        tj <- test_idx[[i]][j]
        x2 <- list(rbind(xtrain[[i]], xtest[[i]][xtest[[i]][,1] <= tj, ]))
        x2_neighbors <- sapply(neighbors$neighbors[[i]],
                               function(n) rbind(xtrain[[n]],xtest[[n]][xtest[[n]][,1] <= tj, ]), simplify=FALSE)
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
                model(x, ...),
              models = c(models[i], models[neighbors$neighbors[[i]]]),
              h = h
            ),
            ndata
          )
        }
        else {
          tmp_res <- forecast(refit(models[[i]], x2[[1]][,-1]), h = h)
        }
        result_i[j+1, , ] <- cbind(max(x2[[1]][,1])+(1:h), tmp_res)
      }
    }
    result <- c(result, list(result_i))
  }

  if (is.null(xtest_idx)){
    result <- lapply(result, "[",1,TRUE,TRUE,drop=FALSE)
    cum_test_errors <- test_errors <- NULL
  }
  else {
    cum_test_errors <- sapply(seq_along(result), function(i) {
      apply(result[[i]], 2, function(r) {
        ind <- intersect(r[,1], xtest[[i]][,1])
        if (length(ind) == 0) return(NA)
        scaling <- if (benchmark == "rw")
          rmse(xtrain[[i]], rw(xtrain[[i]]))
        else
          rmse(cbind(0, residuals(models[[i]])), 0)
        rmsse(xtest[[i]][xtest[[i]][,1] %in% ind, , drop=FALSE], r[r[,1] %in% ind, , drop=FALSE],
              scaling = scaling,
              cumulative = TRUE)
      })
    }, simplify=FALSE)

    test_errors <- sapply(seq_along(result), function(i) {
      apply(result[[i]], 2, function(r) {
        ind <- intersect(r[,1], xtest[[i]][,1])
        if (length(ind) == 0) return(NA)
        scaling <- if (benchmark == "rw")
          rmse(xtrain[[i]], rw(xtrain[[i]]))
        else
          rmse(cbind(0, residuals(models[[i]])), 0)
        rmsse(xtest[[i]][xtest[[i]][,1] %in% ind, , drop=FALSE], r[r[,1] %in% ind, , drop=FALSE],
              scaling = scaling,
              cumulative = FALSE)
      })
    }, simplify=FALSE)
  }


  res <- list(data = x,
              type = type,
              forecasts = result,
              h=h,
              xtest_idx=xtest_idx,
              test_errors = test_errors,
              cum_test_errors = cum_test_errors,
              neighbors = neighbors)

  attr(res,"class") <- "ts_avg"

  return(res)
}

#' plot.ts_avg
#'
#' Plot ts_avg result
#' @param object result of ts_avg function
#' @param plot_type response or error
#' @param h forecast horizon to use. default=1.
#' @return ggplot2 object
#' @example demo/demo.R
#' @import ggplot2
#' @import reshape2
#' @export

plot.ts_avg <- function(object, plot_type=c("response","error"), h=1){

  plot_type <- match.arg(plot_type)
  if (plot_type == "response") {
    x_df <- cbind(avg_type=NA,
                  forecast = FALSE,
                  as.data.frame(do.call(
                    rbind, sapply(seq_along(object$data), function(i)
                      cbind(id = i, fc_idx = NA, object$data[[i]]),
                      simplify=FALSE)
                  )))

    h_plot <- if (object$h == 1) {
      warning(sprintf("parameter h=%d is not used", h))
      1
    }
    else {
      h
    }

    fc_df <- cbind(avg_type=object$type,
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

    plot_df <- rbind(x_df, fc_df)

    colnames(plot_df)[1:5] <- c("avg_type", "forecast", "id", "fc_idx", "index")
    colnames(plot_df)[-(1:5)] <- paste0("V", 1:(ncol(plot_df) - 5))

    plot_res(plot_df, h_plot)
  }
  else if (plot_type == "error"){
    error_df <- as.data.frame(cbind(1:length(object$test_errors), do.call(rbind,object$test_errors)[,h]))
    #error_df <- error_df[order(error_df[,2]),]
    colnames(error_df) <- c("id","rmsse")
    error_df$id <- factor(error_df$id, levels=error_df[order(error_df[,2]),1])

    ggplot2::ggplot(error_df, ggplot2::aes(x=id,y=rmsse, group=1))+
      ggplot2::geom_step()+
      ggplot2::geom_hline(yintercept=1,linetype="dashed")
  }
}
