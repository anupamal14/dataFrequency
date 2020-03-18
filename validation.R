
Acc <- function(err,Y) {
  #------------------------------------------------------
  # Compute and return error metrics:
  # Mean Squared Error (MSE)
  # Mean Absolute Deviation
  # Mean Absolute Percentage Error
  #------------------------------------------------------
  MSE <- round(mean(err^2),2)
  #MAD <- round(mean(abs(err)),2)
  MAPE <- round(mean(abs(err/Y)),4) * 100
  ErrVec <- cbind(MSE, 
                  #MAD, 
                  paste(toString(MAPE),"%"))
  colnames(ErrVec) <- c("MSE", "MAPE") #c("MSE","MAD","MAPE")
  return(ErrVec)
}
