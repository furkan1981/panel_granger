EK <- function(series, NN, dmax = rep(1,NN), lag.max = NULL, ic = c("AIC","BIC"), type = c("constant", "trend")) {
  
  series <- as.matrix(series)
  TT <- nrow(series)/NN             
  nvar <- ncol(series)
  
  ind_lag <- rep(0, NN)
  ind_Wald <- array(0, c(nvar,nvar,NN))
  ind_pvalue <- array(0, c(nvar,nvar,NN))
  
  indix1 <-  matrix(1:(TT*NN),TT,NN) 
  
  for (i in 1:NN) {
    ind_series <- series[indix1[,i],]
    model1 <- TYGC(ind_series, dmax = dmax[i], lag.max, ic , type) 
    ind_lag[i] <-  model1$p
    ind_Wald[,,i] <- model1$Wald  
    ind_pvalue[,,i] <- model1$p_value
  }
  
  Fisher_Stat <- matrix(NA, nvar, nvar)
  p_value <- matrix(NA, nvar, nvar)
  
  indix2 <- matrix(1:nvar,nvar,nvar, byrow = TRUE)
  diag(indix2) <- NA
  indix2 <- t(matrix(t(indix2)[which(!is.na(indix2))], nrow = nvar-1, ncol = nvar))

  for (i in 1:nvar) {
    for (j in indix2[i,]) {
      Fisher_Stat[i,j] <- -2*sum(log(ind_pvalue[i,j,]))
      p_value[i,j] <- pchisq(Fisher_Stat[i,j], df = 2*NN, lower.tail = FALSE)
    }
  }
  
  list("Individual_Wald" = ind_Wald, 
       "Individual_pvalue" = ind_pvalue,
       "Individaul_dmax" = dmax,
       "Individual_LagOrder" = ind_lag,
       "Fisher_Stat" = Fisher_Stat,
       "p_value" = p_value,
       "ic" = ic,
       "type" = type,
       "N" = NN, 
       "T" = TT)
}

