rmsse <- function(y, yhat, benchmark=rw(y), scaling=NULL, cumulative=FALSE){
  scaling_benchmark <- if (is.null(scaling)){
    matrix(sapply(1:nrow(y), function(i) rmse(y[1:i, ,drop=FALSE], benchmark[1:i, ,drop=FALSE])), ncol=ncol(y)-1, nrow=nrow(y),byrow=TRUE)
  }
  else {
    matrix(scaling,ncol=ncol(y)-1,nrow=nrow(y),byrow = TRUE)
  }
  
  scaled_errors <- cbind(y[,1],(y-yhat)[,-1,drop=FALSE]/scaling_benchmark)
  
  if (cumulative) {
    return(sapply(1:nrow(y), function(i) mean(rmse(scaled_errors[1:i, , drop=FALSE],0))))
  }
  else {
    return(mean(rmse(scaled_errors,0)))
  }
}

rmse <- function(y, yhat){
  sqrt(colMeans( (y-yhat)[,-1,drop=FALSE]^2, na.rm=TRUE))
}

rw <- function(y){
  cbind(y[,1],rbind(NA,head(y[,-1,drop=FALSE],-1)))
}

plot_res <- function(plot_df, h){
  if (ncol(plot_df) > 6) {
    plot_df <- reshape(
      plot_df,
      direction = "long",
      varying = 6:ncol(plot_df),
      sep = "",
      idvar = c("forecast", "id", "index"),
      v.names = "value",
      timevar = "comp"
    )
    p <-
      ggplot2::ggplot(data = plot_df, ggplot2::aes(
        x = index,
        y = value,
        linetype = factor(comp)
      ))
  }
  else {
    p <- ggplot2::ggplot(data = plot_df, ggplot2::aes(x = index, y = V1))
  }
  
  p +
    ggplot2::geom_line(ggplot2::aes(color = avg_type)) +
    ggplot2::geom_point(ggplot2::aes(shape = forecast, color=avg_type)) +
    ggplot2::facet_wrap( ~ id, ncol = 4, scales = "free") +
    ggplot2::geom_vline(ggplot2::aes(xintercept = fc_idx), linetype = "dotted") +
    ggplot2::ggtitle(sprintf("Horizon = %d", h))
}