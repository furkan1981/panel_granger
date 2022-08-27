EK_Boot <- function(series, NN, dmax = rep(1,NN), p, type = c("constant", "trend"), 
                    from = 2, to = 1, nboot = 1000) {
  
  series <- as.matrix(series)
  TT <- nrow(series)/NN             
  nvar <- ncol(series)
  
  indix2 <- matrix(1:nvar,nvar,nvar, byrow = TRUE)
  diag(indix2) <- NA
  indix2 <- t(matrix(t(indix2)[which(!is.na(indix2))], nrow = nvar-1, ncol = nvar))
  
  indix1 <-  matrix(1:(TT*NN), TT, NN) 
  
  uhat <- array(NA, c(TT,NN,nvar))
  Ahat <- NULL
  
  if (type == "constant") l = 1  else if (type == "trend") l = 2
  
  m <- max(p) + max(dmax) 
  Ahat <- array(NA, c(nvar, m*nvar+l, NN))
  
  for (i in 1:NN) {
    ind_series <- series[indix1[,i],]
    if (type == "constant") {
      model_var <- VAR(ind_series, p[i] + dmax[i], type = "const")
    } else if (type == "trend") {
      model_var <- VAR(ind_series, p[i] + dmax[i], type = "both")
    }
  
    A <- matrix(1, nrow = nvar, ncol = nvar*(p[i] + dmax[i])+l)
    
    for (j in 1:p[i]) A[to,from*j] = 0
    
    model_var_rest <- restrict(model_var, method = "manual", resmat = A)
    
    u_ind <- residuals(model_var_rest)
    for (k in 1:nvar) uhat[(p[i]+dmax[i]+1):TT,i,k] <- u_ind[,k] 
    
    Ahat[1:nvar,1:((p[i]+dmax[i])*nvar+l),i] <- Bcoef(model_var_rest)
  }
  
  uhat <- uhat[-(1:(max(p)+max(dmax))),,]
  
  uhattilda <- array(NA, c(nrow(uhat),NN,nvar))
  for (k in 1:nvar) uhattilda[,,k] <- demean(uhat[,,k]) 

  # starting bootstrap
  Fisher_Stat_Boot <- NULL
  
  for (s in 1:nboot) {
    star <- sample(1:nrow(uhattilda), size = TT, replace = TRUE) 
    ustar <- uhattilda[star,,]
    
    ystar <- matrix(NA, ncol = nvar, nrow = NN*TT)
    
    for (i in 1:NN) {
      Ahat1 <- Ahat[,1:((p[i]+dmax[i])*nvar),i]
      Ahat2 <- Ahat[,(ncol(Ahat1)+1):(ncol(Ahat1)+l),i]
      Ahat3 <- cbind(Ahat2, Ahat1)
      ustarnew <- ustar[,i,]
      
      if (type == "constant") {
        ystar[indix1[,i],] <- VAR.sim(B = Ahat3, n = TT, lag = p[i]+dmax[i], innov = ustarnew, include = "const")  
      } else if (type == "trend") {
        ystar[indix1[,i],] <- VAR.sim(B = Ahat3, n = TT, lag = p[i]+dmax[i], innov = ustarnew, include = "both")  
      }
    }
    
    ind_Wald <- array(0, c(nvar,nvar,NN))
    ind_pvalue <- array(0, c(nvar,nvar,NN))
    
    for (i in 1:NN) {
      ind_series <- ystar[indix1[,i],]
      model1 <- TYGC_Boot(ind_series, dmax[i], lag = p[i], type) 
      ind_Wald[,,i] <- model1$Wald  
      ind_pvalue[,,i] <- model1$p_value
    }
    
    Fisher_Stat_Boot[s] <- -2*sum(log(ind_pvalue[to,from,]))
  }
  
  CV_Boot <- quantile(Fisher_Stat_Boot, c(.99,.95,.90))
  
  list("from_variable" = from,
       "to_variable" = to,
       "Bootstrap_Critical_Values" = CV_Boot)
  }
  

    
    

    

   

  
  
  
  


