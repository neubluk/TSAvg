#' ts_avg.cv
#'
#' Perform TS Avg Forecasting TS Cross-Validation for specific Avg. Types
#'
#' The forecast horizon for TSCV and resulting best models is fixed to h=1.
#'
#' @param x list of matrices. each matrix corresponds to one time series whereby the first column always corresponds to its index.
#' @param model a generic model function. the resulting model object must have forecast, refit, fitted, and residuals functions (as demonstrated in my_ets).
#' @param min_ts smallest time index to do TSCV, i.e. the first fold contains all time series cut off to 1:min_ts.
#' @param k_grid cv grid for k.
#' @param type averaging type. subset of  "simple_avg", "simple_avg_n", "dist_avg_n", "dist_avg", "global_avg", "perf_avg", "perf_avg_r". (default is all).
#' @param xtest_idx index to split the every time series in train and test set. default = NULL.
#' @param benchmark benchmark forecast type for computing the errors. either random walk forecast "rw" or model forecast "model".
#' @param verbose if TRUE, then progress is printed.
#' @param ... additional parameters for model function
#' @return object of class ts_avg.cv
#' @example demo/demo.R
#' @seealso [plot.ts_avg.cv]
#' @export
ts_avg.cv <- function(x, model, min_ts, k_grid, type=c("simple_avg","simple_avg_n", "dist_avg_n", "dist_avg","global_avg","perf_avg","perf_avg_r"), xtest_idx = NULL, benchmark=c("rw","model"), verbose=FALSE, ...){

  timesteps <- sort(unique(unlist(lapply(x,"[",TRUE,1))))

  timesteps <- timesteps[timesteps>=min_ts]

  if (!is.null(xtest_idx)){
    timesteps <- timesteps[timesteps < xtest_idx]
  }

  result <- rmsses <- array(dim = c(length(x),length(timesteps),ncol(x[[1]]),length(type),length(k_grid)))

  cv_scores <- cv_ses <- array(dim=c(length(x),length(type),length(k_grid)),
                               dimnames=list(id=seq_along(x),type=type,k=k_grid))

  for (i in seq_along(type)) {
    for (j in seq_along(k_grid)) {
      for (l in seq_along(timesteps)) {
        x_train <- lapply(x, function(xi)
          xi[xi[, 1] < timesteps[l], , drop=FALSE])

        x_test <- lapply(x, function(xi)
          xi[xi[, 1] >= timesteps[l] - l + 1 & xi[,1] <= timesteps[l], , drop=FALSE])

        x_train_lengths <- unlist(lapply(x_train,nrow))

        result[x_train_lengths <= 1, l, , i, j] <- NA
        if (verbose) message(sprintf("Performing avg. type %s with k=%d up to index %d",type[i],k_grid[j],timesteps[l]))
        result_fc <- ts_avg(x_train[x_train_lengths>1], model, k_grid[j], h=1, type[i], xtest_idx = NULL, ...)

        for (m in which(x_train_lengths>1)){
          if (max(x[[m]][, 1]) >= timesteps[l]) {
            result[m, l, , i, j] <- result_fc$forecasts[[which(which(x_train_lengths > 1) == m)]]

            rmsses[m, l, , i, j] <-
              rmsse(x_test[[m]], result[m, (l-nrow(x_test[[m]])+1):l, , i, j])
          }
        }
      }
      cv_scores[, i, j] <- apply(rmsses[, , , i, j, drop=FALSE], 1, mean, na.rm = TRUE)
      cv_ses[, i, j] <- apply(rmsses[, , , i, j, drop=FALSE], 1, function(z) sd(z,na.rm = TRUE)/sqrt(sum(!is.na(z))))
    }
  }

  best_idx <- apply(cv_scores,1:2,function(z){
    if (all(is.na(z))) return(0)
    which.min(z)
  })
  # best_idx_se <- apply(cv_scores <= array((cv_scores+cv_ses)[best_idx],dim=c(length(x),length(type),length(k_grid))),1:2,function(z){
  #   if (all(is.na(z)) | all(!(z))) return(NA)
  #   min(which(z))
  # })
  best_idx_se <- apply(expand.grid(dimnames(best_idx)),
                       1,
                       function(z) min(which( cv_scores[z[1],z[2],] <= (cv_scores+cv_ses)[z[1],z[2],best_idx[z[1],z[2]]])))

  best_idx_se[is.infinite(best_idx_se)] <- best_idx[is.infinite(best_idx_se)]

  best_idx_final <- (best_idx + best_idx_se - abs(best_idx-best_idx_se))/2

  colnames(best_idx_final) <- type

  best_k <- apply(best_idx_final, 1:2, function(z){
    if (z > 0) return(k_grid[z])
    return(0)
  })

  best_result <- lapply(type,function(t){
    ts_avg(x, model, best_k[,t], h=1, type=t, xtest_idx = xtest_idx,benchmark=benchmark)
  })

  names(best_result) <- type

  best_models <- apply(simplify2array(lapply(best_result,'[[','test_errors')),1,which.min)

  final_best_models <- list()
  for (i in seq_along(x)){
    final_best_models <- c(final_best_models, best_result[best_models[i]])
  }

  res <- list(best_result=best_result,
              best_k=best_k,
              final_models = final_best_models,
              cv_results=list(scores=cv_scores,
                              se=cv_ses))
  attr(res,"class") <- "ts_avg.cv"

  return(res)
}

#' plot.ts_avg.cv
#'
#' Plot ts_avg.cv result
#' @param object result of ts_avg.cv function
#' @param plot_type cv results, or best forecasts, or error
#' @param type select only certain avg. types for plotting
#' @return ggplot2 object
#' @example demo/demo.R
#' @import ggplot2
#' @import reshape2
#' @export
#'
plot.ts_avg.cv <- function(object, plot_type=c("cv","best","error"), type=NULL){
  # todo
  # plot errors of best models

  plot_type <- match.arg(plot_type)
  if (plot_type == "cv"){
    plot_df <- cbind(reshape2::melt(res_mv$cv_results$scores, value.name = "cverror"),
                     se=as.numeric(res_mv$cv_results$se))

    plot_df <- merge(plot_df,
                     reshape2::melt(object$best_k, value.name = "best_k"))

    plot_df$best_k[plot_df$best_k==0] <- NA

    if (!is.null(type)){
      plot_df <- plot_df[plot_df$type %in% type, ]
    }

    ggplot2::ggplot(plot_df, ggplot2::aes(x=k, color=type))+
      ggplot2::geom_line(ggplot2::aes(y=cverror))+
      ggplot2::geom_errorbar(ggplot2::aes(ymin=cverror-se,ymax=cverror+se))+
      ggplot2::geom_vline(ggplot2::aes(xintercept=best_k, color=type), linetype="dashed")+
      ggplot2::facet_wrap(~id, ncol=4, scales="free")
  }
  else if (plot_type == "best"){
    plot_df <- sapply(seq_along(object$final_models), function(i){
      fc_df <- as.data.frame(object$final_models[[i]]$forecasts[[i]][,1,])
      rbind(cbind(avg_type=NA,forecast=FALSE,id=i,fc_idx=min(fc_df[,1])-1,as.data.frame(object$final_models[[i]]$data[[i]])),
            cbind(avg_type=object$final_models[[i]]$type,forecast=TRUE,id=i,fc_idx=min(fc_df[,1])-1,fc_df))
    }, simplify=FALSE)

    plot_df <- do.call(rbind, plot_df)
    colnames(plot_df)[5] <- "index"
    colnames(plot_df)[-(1:5)] <- paste0("V", 1:(ncol(plot_df)-5))

    plot_res(plot_df,h=1)
  }
  else if (plot_type == "error"){
    plot_df <- stack(lapply(res_mv$best_result, function(b) unlist(b$cum_test_errors)))

    ggplot2::ggplot(plot_df, ggplot2::aes(x=ind, y=values))+
      ggplot2::geom_boxplot()+
      ggplot2::geom_hline(yintercept=1)
  }
}
