source("inst/my_ets.R")

x <- replicate(5,{
  len <- rpois(1,20)
  cbind(seq(rpois(1,5),length.out=len),matrix(rnorm(2*len),ncol = 2))
  },simplify=FALSE)


ts_avg_res <- ts_avg(x,my_ets,h=2,k=2,type="simple_avg", xtest_idx = 20)
plot(ts_avg_res)

ts_avg_res_2 <- ts_avg(lapply(x,"[",TRUE,1:2),my_ets,h=2,k=1,type="simple_avg", xtest_idx = 15)
plot(ts_avg_res_2,h=1)

# cv examples
res_mv <- ts_avg.cv(x,my_ets,min_ts=10,k_grid=2:5,xtest_idx = 20,verbose=TRUE)

plot(res_mv, "cv")
plot(res_mv, "error")
plot(res_mv, "best")
