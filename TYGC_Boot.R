TYGC_Boot <- function(series, dmax = 1, lag = NULL, type = c("constant","trend")) { 
  series <- as.matrix(series)
  nvar <- ncol(series)        # Number of variables
  n <- nrow(series)           # Number of observations
  
  trend <- c(1:n)
  
  seriesdep <- series[-(1:(lag + dmax)),]
  serieslag <- embed(series, lag + dmax + 1)[,-(1:nvar)]
  
  test_stat <- matrix(0,nvar,nvar)
  p_value <- matrix(0,nvar,nvar)
  
  R <- matrix(2:(nvar*lag + 1),nvar,lag)
  
  for (i in 1:nvar) {
    if (type == "constant") {
      result <- lm(seriesdep[,i] ~ serieslag)
    } else if (type == "trend") {
      result <- lm(seriesdep[,i] ~ serieslag + trend[-(1:(lag+dmax))])
    } 
    
    for (j in 1:nvar) {
      wald_result <- wald.test(b = coef(result), Sigma = vcov(result), Terms = R[j,])
      test_stat[i,j] <- wald_result$result$chi2[["chi2"]]
      p_value[i,j] <- wald_result$result$chi2[["P"]]
    }
  }
  list("Wald" = test_stat, 
       "p_value" = p_value,
       "dmax" = dmax,
       "type" = type)
}  




