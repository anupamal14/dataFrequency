tbatsHighFreq <- function(x,Start,seas,steps,n){
  N <- length(x)
  # Multiple Seasonality series
  xMSTS <- msts(x[1:n],start=Start,seasonal.periods=seas)
  #Fit TBATS model
  fit <- tbats(xMSTS)
  # Model period values
  xhat<- fit$fitted.values
  yhat <- as.vector(predict(fit,h=steps)$mean)
  err <- list(model = (x[1:n] - xhat), hold = (x[(n+1):(n+steps)] - yhat))
  return(err)
}