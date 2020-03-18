WKestimatorPA <- function(ARMAmodel, data){
  #------------------------------------------------------
  # Wiener Kolmogorov estimator for 
  c0 <- ARMAmodel$coef[1]
  d <- ARMAmodel$coef[3]
  a <- -ARMAmodel$coef[2] # MA part negated to be consistent with +signs in the filter defn.
  b <- -ARMAmodel$coef[4] # MA part negated to be consistent with +signs in the filter defn.
  
 m <- mean(data)
 
  centeredData <- data - m
  ywhitened <- stats::filter(centeredData, c(0,0,0,b), method = "recursive")
  ywhitened <- stats::filter(ywhitened, a, method = "recursive")
  
  ywhitened <- stats::filter(ywhitened, c(1,-c0), method = "convolution", 
                             sides = 1)
  ywhitened <- stats::filter(ywhitened, c(1,0,0,0,-d), method = "convolution", 
                             sides = 1)

  # Now, compute the error sequence, delayed by 4
  m3 <- 1 - a + c0 + (-a*c0 + c0^2) + (-a*c0^2 + c0^3)
  m2 <- 1 - a + c0 + (-a*c0 + c0^2)
  m1 <- 1 - a + c0
  m0 <- 1
  yerr <- stats::filter(ywhitened, c(m0,m1,m2,m3), method = "convolution", 
                        sides = 1)
  
  # Delay by one more, because each convolution above (there are three) already yields
  # a delay of one, and the fourth here ensures that the error at t is the prediction error
  # for the sum y(t)+y(t-1)+y(t-2)+y(t-3)
  yerr <- stats::filter(yerr, 1, method = "convolution", 
                        sides = 1)
  
  # The series is delayed by 4. Also next 4 are NAs. So, fill up these with 
  # prediction errors arising from predicting the value to be the mean
  for (k in (1:8)){
    yerr[k] <- sum(data[max(1,k-3):k]) - min(k,4)*m
  }
  
  # Return every 4th sample
  yerr_sub4 <- yerr[seq(4,length(yerr), 4)]
  
  return(yerr_sub4)
}

WKestimatorAP <- function(ARMAmodel, data) {
  c <- ARMAmodel$coef[1]
  d <- ARMAmodel$coef[3]
  a <- -ARMAmodel$coef[2] # MA part negated to be consistent with +signs in the filter defn.
  b <- -ARMAmodel$coef[4] # MA part negated to be consistent with +signs in the filter defn.
  
  
  denom <- (-a*((1+c)^2)*(1+4*c^2 + c^4) + (1 + a^2)*(c + 2*c^2 + 4*c^3 + 2*c^4 + c^5))
  num <- (-a*((1 + c)^2)*(6 + 8*c^2 + 6*c^4) + (1 + a^2)*
            (4 + 6*c + 8*c^2 + 8*c^3 + 8*c^4 + 6*c^5 + 4*c^6))
  alpha <- -(num - sqrt(num^2 - 4*(denom^2)))/(2*denom)
  # Keeping gamma always positive, hence abs(denom).
  gamma <- sqrt(2)*abs(denom)/sqrt(num - sqrt(num^2 - 4*(denom^2)))
  
  centeredData <- data - mean(data)
  ywhitened <- stats::filter(centeredData, alpha, method = "recursive")
  ywhitened <- stats::filter(ywhitened, b, method = "recursive")
  
  ywhitened <- stats::filter(ywhitened, c(1,-d), method = "convolution", 
                             sides = 1)
  ywhitened <- stats::filter(ywhitened, c(1,-(c^4)), method = "convolution", 
                             sides = 1)
  
  # Convolution gives two delays. Undo one of them by advancing by one sample.
  yerr <- ywhitened[2:length(ywhitened)]
  yerr[1] <- data[1] - mean(data)
  
  return(yerr)  
}
