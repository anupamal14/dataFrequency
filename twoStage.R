multSeasonality <- function(data,BlockFreq, NrDaysPerWk, NrHrsInADay) {
  
  #BlockFreq   <- 12
  #NrHrsInADay <- 24
  DayFreq     <- NrHrsInADay*BlockFreq
  #NrHrsInADay      <- 7
   
  if (NrDaysPerWk == 7) {
    dayNames <- c("sun","mon","tue","wed","thu","fri","sat")
  } else {
    dayNames <- c("mon","tue","wed","thu","fri")
  }  
  n <- nrow(data)
  
  #------------------------------------------------------
  # Block seasonality
  #------------------------------------------------------
  d_hour <- colMeans(matrix(data[,1], nrow = BlockFreq))
  # Replicate this 12 times and reshape into a vector so that 
  # it is the same length as d
  d_hour_rep <- as.vector(t(matrix(t(rep(t(d_hour), BlockFreq)), ncol = BlockFreq)))
  # Get block seasonal index
  b.ratio <- data[,1]/d_hour_rep
  b.day <- rowSums(data[,10:16]*matrix(rep(c(1:NrDaysPerWk),n),nrow=n,byrow=TRUE))
  tmp <- rep(c(1:DayFreq),ceiling(n/DayFreq))
  b.block <- tmp[1:n]
  
  a <- aggregate(b.ratio,list(day = b.day, block = b.block),mean)
  block <- matrix(a[,3],nrow = DayFreq, ncol = NrDaysPerWk, byrow=TRUE)
  colnames(block) <- dayNames
  
  #--------------------------------------------------------------------------
  # Hourly seasonality indices
  #--------------------------------------------------------------------------
  # First get daily averages and then divide hourly averages by daily averages
  d_day <- colMeans(matrix(data[,1], nrow = DayFreq))
  # Replicate by 24 so that we have the same value for each hour in a day
  d_day_rep <- as.vector(t(matrix(t(rep(t(d_day), NrHrsInADay)), ncol = NrHrsInADay)))
  # Get hourly seasonal index
  b2.ratio <- d_hour/d_day_rep
  tmp <- rep(1:NrHrsInADay,ceiling(length(d_hour)/NrHrsInADay))#as.vector(t(matrix(rep(1:NrHrsInADay,ceiling(length(d_hour)/NrHrsInADay)), nrow=NrHrsInADay)))
  b2.hour <- tmp[1:length(d_hour)]
  b2.day <- b.day[seq(1,length(b.day),DayFreq/NrHrsInADay)]
  a2 <- aggregate(b2.ratio,list(day = b2.day, hour = b2.hour),mean)
  hourly <- matrix(a2[,3],nrow = NrHrsInADay, byrow = TRUE)
  colnames(hourly) <- dayNames
  
  
  # Initialize list to return the values
  SI <- list()
  SI$block  <- block
  SI$hourly <- hourly
  #SI$daily  <- daily
  return(SI)
  
}


# Remove the variables (covariates) which are not significant
bestModel <- function(d,y,fit) {
  s <- summary(fit)
  # Coefficients which are not significant at 10% p-value
  idx <- which(s$coefficients[,4] > 0.1)
  if (length(idx) == 0) {
    bestFit <- fit
  } else {
    idx2 <- which(s$coefficients[idx,4] == max(s$coefficients[idx,4]))
    str11 <- toString(c(rownames(s$coefficients)[c(-1,-idx[idx2])]))
    x <- gsub(","," +",str11)
    str2 <- paste("lm(",y," ~ ",x,", data = d)")
    bestFit <- eval(parse(text=str2))
    s <- summary(bestFit)
    idx <- which(s$coefficients[,4] > 0.1)
    idx2 <- which(s$coefficients[idx,4] == max(s$coefficients[idx,4]))
    if (length(idx) != 0) {
      idxToRemove <- which(colnames(d) == names(idx[idx2])) 
      bestFit <- bestModel(d[,-idxToRemove],y,bestFit) 
    }
  }
  return(bestFit)
}


twoStageOneStep <- function(d, n, seas, NrDaysPerWk, NrHrsPerDay){
  d <- cbind(d,d[,c(19:22)]^2)
  colnames(d)[19:26] <- c("tmax","tmin","bmax","bmin","tmaxSq","tminSq","bmaxSq","bminSq")
  colnames(d)[9] <- "DNo"
  nrSampsSum <- 96
  lambda <- 4  
  
  
  dailyData <- d[seq(1,nrow(d),nrSampsSum),]
  dailyData$load <- colSums(matrix(d$load, nrSampsSum, length(d$load)/nrSampsSum))/nrSampsSum
  
  DL <- dailyData[1:(n/nrSampsSum),]
  
  # Find the start day index of the test period
  testPeriodDayStartInd <- d$dayweek[n+1]
  
  
  fitD2 <- lm(load ~ DNo + tmax + tmin + bmax + bmin +
                tmaxSq + tminSq + bmaxSq + bminSq + weekday, data = DL)
  
  # Include weekday and saturday
  fitF2 <- lm(load ~ DNo + tmax + tmin + bmax + bmin +
                tmaxSq + tminSq + bmaxSq + bminSq + weekday + sat, data = DL)
  
  # Model H with square temp variables
  fitH2 <- lm(load ~ DNo + tmax + tmin + bmax + bmin + 
                tmaxSq + tminSq + bmaxSq + bminSq +
                mon + tue + wed + thu + fri + sat, data = DL)
  
  a1 <- anova(fitD2,fitF2)
  pval <- a1$`Pr(>F)`
  if (pval[2] > .05)  {bmodel <- fitD2 } else {bmodel <- fitF2}
  a2 <- anova(bmodel,fitH2)
  pval <- a2$`Pr(>F)`
  if(pval[2] < .05) { bmodel <- fitH2 }
  
  bestFit <- bestModel(DL,"load",bmodel)
  
  
  xMSTS <- msts(bestFit$residuals,start=8+248/365,seasonal.periods=c(7,365))
  
  
  fit <- tbats(xMSTS)
  
  
  xhat2 <- bestFit$fitted.values + fit$fitted.values
 
 
  # Predict and aggregate
  #--------------------------------------
  xhatA <- as.vector(matrix(rep(xhat2,nrSampsSum), nrow=nrSampsSum, byrow=TRUE))
  # Just regression
  x2A <- as.vector(matrix(rep(bestFit$fitted.values,nrSampsSum), nrow=nrSampsSum, byrow=TRUE))
  #-----------------------------------------
  
  #-----------------------------------------
  # Average and predict
  #-----------------------------------------
  # Replicate to get 24 values for each day
  xhat <- as.vector(matrix(rep(xhat2, nrSampsSum/lambda), nrow=nrSampsSum/lambda, byrow=TRUE))
  # Just regression
  x2 <- as.vector(matrix(rep(bestFit$fitted.values, nrSampsSum/lambda), nrow=nrSampsSum/lambda, byrow=TRUE))
  #-----------------------------------------
  
  rowsIndices2 <- (n+1):(nrow(d))#(n+lambda)
  y <- d$load[rowsIndices2]
  yHr <- colMeans(matrix(y,nrow = lambda))
  xHr <- colMeans(matrix(d$load[1:n],nrow = lambda))
  rowsIndices <- n/nrSampsSum + (1:(length(yHr)/24))
  
  colsForNewData <- matrix(NA,nrow=1,ncol=(length(rownames(summary(bestFit)$coefficients))-1))
  for (i2 in c(2:length( rownames(summary(bestFit)$coefficients)))) {
    colsForNewData[i2-1] <-  which(colnames(dailyData) == rownames(summary(bestFit)$coefficients)[i2])
  }
  newRegData <- dailyData[rowsIndices,colsForNewData]
  y1<- predict(bestFit,newdata = newRegData)
  yhat2 <- y1 + as.vector(predict(fit,h=1)$mean)
  
  #----------------------------------------
  # Predict and Aggregate
  yhatA <- as.vector(matrix(rep(yhat2,nrSampsSum),nrow=nrSampsSum,byrow=TRUE))
  # Just regression
  y2A <- as.vector(matrix(rep(y1,nrSampsSum),nrow=nrSampsSum,byrow=TRUE))
  #-----------------------------------------
  
  #-----------------------------------------
  # Average and predict
  #-----------------------------------------
  # Replicate to get 24 values for each day
  yhat <- as.vector(matrix(rep(yhat2,nrSampsSum/lambda), nrow=nrSampsSum/lambda, byrow=TRUE))
  # Just regression
  y2 <- as.vector(matrix(rep(y1,nrSampsSum/lambda), nrow=nrSampsSum/lambda, byrow=TRUE))
  
  
  # Seasonality correction
  #-----------------------
  # Get the seasonal corrections for block and hour because we are not using them in the TBATS call
  # and apply to the predicted values
  si <- multSeasonality(d,lambda, NrDaysPerWk, NrHrsPerDay)
  
  # Hourly seasonal corrections
  SImult1 <- si$hourly[,c(testPeriodDayStartInd:7, 1:(testPeriodDayStartInd-1))]
  SImult_vec1 <- as.vector(SImult1)
  repNr <- ceiling(length(yhat)/length(SImult_vec1))
  SImultHr_vec <- rep(SImult_vec1,times = repNr)
  
  
  # Apply the seasonality corrections
  yhatHr <- yhat*SImultHr_vec[1:length(yhat)]
  # Apply seasonality correcction for plain regression model. 
  # Assuming daily corrections are taken care of by dummy variables for weekdays etc
  y2Hr <- y2*SImultHr_vec[1:length(yhat)]
  
  # Hourly seasonal corrections
  SImult1 <- si$hourly[,c(d$dayweek[1]:7, 1:(d$dayweek[1]-1))]
  SImult_vec1 <- as.vector(SImult1)
  repNr <- ceiling(length(xhat)/length(SImult_vec1))
  SImultHr_vec <- rep(SImult_vec1,times = repNr)
  
  
  # Apply the seasonality corrections
  xhatHr <- xhat*SImultHr_vec[1:length(xhat)]
  # Apply seasonality correcction for plain regression model. 
  # Assuming daily corrections are taken care of by dummy variables for weekdays etc
  x2Hr <- x2*SImultHr_vec[1:length(xhat)]
  #-------------
  
  # Now need both block and hourly corrections
  SImult1 <- si$block[,c(testPeriodDayStartInd:7, 1:(testPeriodDayStartInd-1))]
  SImult_vec1 <- as.vector(SImult1)
  repNr <- ceiling(length(yhatA)/length(SImult_vec1))
  SImultBlk_vec <- rep(SImult_vec1,times = repNr)
  
  SImult1 <- si$hourly[,c(testPeriodDayStartInd:7, 1:(testPeriodDayStartInd-1))]
  SImult_vec1 <- as.vector(SImult1)
  SImulttmp <- as.vector(t(matrix(rep(SImult_vec1, lambda), ncol=lambda)))
  SImultHr_vec <- rep(SImulttmp, repNr)
  
  # Apply the seasonality corrections
  yhatAHr <- colMeans(matrix(yhatA*SImultBlk_vec[1:length(yhatA)]*SImultHr_vec[1:length(yhatA)],nrow=lambda))
  # Apply seasonality correcction for plain regression model
  y2AHr <- colMeans(matrix(y2A*SImultBlk_vec[1:length(yhatA)]*SImultHr_vec[1:length(yhatA)],nrow=lambda))
  
  # Now need both block and hourly corrections
  SImult1 <- si$block[,c(d$dayweek[1]:7, 1:(d$dayweek[1]-1))]
  SImult_vec1 <- as.vector(SImult1)
  repNr <- ceiling(length(xhatA)/length(SImult_vec1))
  SImultBlk_vec <- rep(SImult_vec1,times = repNr)
  
  SImult1 <- si$hourly[,c(d$dayweek[1]:7, 1:(d$dayweek[1]-1))]
  SImult_vec1 <- as.vector(SImult1)
  SImulttmp <- as.vector(t(matrix(rep(SImult_vec1, lambda), ncol=lambda)))
  SImultHr_vec <- rep(SImulttmp, repNr)
  
  
  # Apply the seasonality corrections
  xhatAHr <- colMeans(matrix(xhatA*SImultBlk_vec[1:length(xhatA)]*SImultHr_vec[1:length(xhatA)],nrow=lambda))
  # Apply seasonality correcction for plain regression model
  x2AHr <- colMeans(matrix(x2A*SImultBlk_vec[1:length(xhatA)]*SImultHr_vec[1:length(xhatA)],nrow=lambda))
  #---------------
  
 
  errPA_H <- rbind((yHr - yhatAHr),(yHr - y2AHr))
  errPA_M <- rbind((xHr - xhatAHr),(xHr - x2AHr))
  idx <- which.min(rowMeans(errPA_H^2))
  
  errAP_H <- rbind((yHr - yhatHr),(yHr - y2Hr))
  errAP_M <- rbind((xHr - xhatHr),(xHr - x2Hr))
  idx2 <- which.min(rowMeans(errAP_H^2))
  
  
  retErr <- list(modelPA = errPA_M[idx,], holdPA = errPA_H[idx,],
                 modelAP = errAP_M[idx2,], holdAP = errAP_H[idx2,])
  return(retErr)
}