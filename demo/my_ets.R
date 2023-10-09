require("forecast")

my_ets <- function(y, ...){
  if (is.matrix(y) && ncol(y)>1) {
    ets_fit <- apply(y, 2, forecast::ets, ...)
    attr(ets_fit, "dim") <- ncol(y)
  } else {
    ets_fit <- list(forecast::ets(y, ...))
    attr(ets_fit, "dim") <- 1
  }
  attr(ets_fit, "class") <- "my_ets"
  attr(ets_fit, "failed") <- FALSE
  return(ets_fit)
}

print.my_ets <- function(object){
  lapply(object, function(o) {
    tmp <- o
    attr(tmp, "class") <- "ets"
    print(tmp)
  })
}

forecast <- function(object, h) UseMethod("forecast")

forecast.my_ets <- function(object, h=1){
  if (attr(object,"failed")) return(NA)
  return(do.call(cbind,lapply(object, function(o){
    tmp <- o
    attr(tmp,"class") <- "ets"
    return(matrix(forecast::forecast(tmp, h=h)$mean,nrow=h))
  })))
}

refit <- function(object, ts, ...) UseMethod("refit")

refit.my_ets <- function(object, ynew, ...){
  ynew <- as.matrix(ynew)
  stopifnot(length(object) == ncol(ynew))
  tmp_refit <- sapply(seq_along(object), function(i) {
    tmp <- object[[i]]
    attr(tmp, "class") <- "ets"

    res <- tryCatch({
      forecast::ets(ynew[,i], model = tmp, use.initial.values = FALSE, ...)
    },
    error=function(e){
      return(NULL)
    })
  }, simplify=FALSE)

  attr(tmp_refit, "class") <- "my_ets"
  attr(tmp_refit, "failed") <- FALSE

  if (any(sapply(tmp_refit,is.null))){
    attr(tmp_refit, "failed") <- TRUE
  }

  return(tmp_refit)
}

fitted.my_ets <- function(object,h=1){
  if (attr(object,"failed")) return(NA)
  return(do.call(cbind, lapply(object, function(o) {
    tmp <- o
    attr(tmp, "class") <- "ets"
    return(matrix(fitted(tmp, h = h),ncol=1))
  })))

}

residuals.my_ets <- function(object,...){
  if (attr(object,"failed")) return(NA)
  return(do.call(cbind, lapply(object, function(o) {
    tmp <- o
    attr(tmp, "class") <- "ets"
    return(matrix(residuals(tmp, type="response"),ncol=1))
  })))
}
