# Toda & Yamamoto (1995) Granger Causality 
# written by Furkan Emirmahmutoglu
# emirfurkan@gmail.com; f.emirmahmutoglu@hbv.edu.tr

TYGC <- function(series, dmax = 1, lag.max = NULL, ic = c("AIC", "BIC"), type = c("constant","trend")) { 
  series <- as.matrix(series)
  nvar <- ncol(series)        # Number of variables
  n <- nrow(series)           # Number of observations
  
  trend <- c(1:n)
  
  if (type == "constant") {
    a <- VARselect(series, lag.max)  
  } else if (type == "trend") {
    a <- VARselect(series, lag.max, type = "both")  
  }
  
  if (ic == "AIC") {
    opt_lag <- a$selection[["AIC(n)"]]
  } else if (ic == "BIC") {
    opt_lag <- a$selection[["SC(n)"]]
  }
  
  seriesdep <- series[-(1:(opt_lag + dmax)),]
  serieslag <- embed(series, opt_lag + dmax + 1)[,-(1:nvar)]
  
  test_stat <- matrix(0,nvar,nvar)
  p_value <- matrix(0,nvar,nvar)
  
  R <- matrix(2:(nvar*opt_lag+1),nvar,opt_lag)
  
  # from j to i 
  for (i in 1:nvar) {
    if (type == "constant") {
      result <- lm(seriesdep[,i] ~ serieslag)
    } else if (type == "trend") {
      result <- lm(seriesdep[,i] ~ serieslag + trend[-(1:(opt_lag+dmax))])
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
       "p" = opt_lag,
       "ic" = ic,
       "type" = type)
}  




