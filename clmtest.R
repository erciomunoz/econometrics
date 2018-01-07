clmtest = function(x) 
  {
############## clmtest() ############################################
# For a reference see:
# Baltagi (2013), Econometric Analysis of Panel Data, 5th edition, pp. 72.
# Written by Ercio Munoz, Graduate Center CUNY
# January 2, 2018

  if (x$args$model == "within") {
    x = update(x, model = "random")
    } else if (x$args$model == "pooling") {
      stop("This test uses GLS one-way residuals, pooled OLS not allowed")
    } 
  
  pdim = pdim(x)
  N = pdim$nT$n
  T = pdim$nT$T
#  balanced = pdim$balanced # true or false
  res = resid(x)
  effect = x$args$effect
  
# Pre-allocating matrices
  I_N    = diag(N)
  Jbar_N = matrix(1/N,N,N)
  I_T    = diag(T)
  Jbar_T = matrix(1/T,T,T)
  E_N    = I_N-Jbar_N
  E_T    = I_T-Jbar_T
  e      = matrix(res,N*T,1)
  et     = matrix(res,1,N*T)
  
  ### START calculation of parts of test statistic ##
#  sigma_2_2   = ( matrix(res,1,N*T) %*% kronecker(matrix(1/N,N,N),diag(T)) %*% matrix(res,N*T,1) ) / T 
  sigma_2_2   = ( crossprod(e , tcrossprod(kronecker(Jbar_N,I_T),et) ) ) / T 
  
#  sigma_2_eta1 = (matrix(res,1,N*T) %*% kronecker(diag(N)-matrix(1/N,N,N),diag(T)) %*% matrix(res,N*T,1)) / (T*(N-1)) 
  sigma_2_eta1 = ( crossprod(e,tcrossprod(kronecker(E_N,I_T),et)) ) / (T*(N-1)) 
  
#  sigma_2_1 = (matrix(res,1,N*T) %*% kronecker(diag(N),matrix(1/T,T,T)) %*% matrix(res,N*T,1)) / N 
  sigma_2_1 = ( crossprod(e,tcrossprod(kronecker(I_N,Jbar_T),et)) ) / N 
  
#  sigma_2_eta2 = (matrix(res,1,N*T) %*% kronecker(diag(N),diag(T)-matrix(1/T,T,T)) %*% matrix(res,N*T,1)) / (N*(T-1)) 
  sigma_2_eta2 = ( crossprod(e , tcrossprod(kronecker(I_N,E_T) , et)) ) / (N*(T-1)) 
  
#  D11 = (matrix(res,1,N*T) %*% kronecker(matrix(1/N,N,N),matrix(1/T,T,T)) %*% matrix(res,N*T,1)) / sigma_2_2
  D11 = ( crossprod(e , tcrossprod(kronecker(Jbar_N,Jbar_T) , et)) ) / sigma_2_2
  
#  D12 = (matrix(res,1,N*T) %*% kronecker(diag(N)-matrix(1/N,N,N),matrix(1/T,T,T)) %*% matrix(res,N*T,1)) / (sigma_2_eta1*(N-1))
  D12 = ( crossprod(e , tcrossprod(kronecker(E_N,Jbar_T) , et)) ) / (sigma_2_eta1*(N-1)) 
   
  D1  = (T/2) * ((1/sigma_2_2)*(D11 - 1) + (N-1) * (1/sigma_2_eta1) * (D12-1) )
  
#  D21 = (matrix(res,1,N*T) %*% kronecker(matrix(1/N,N,N),matrix(1/T,T,T)) %*% matrix(res,N*T,1)) / sigma_2_1
  D21 = ( crossprod(e , tcrossprod(kronecker(Jbar_N,Jbar_T) , et)) ) / sigma_2_1
  
#  D22 = (matrix(res,1,N*T) %*% kronecker(matrix(1/N,N,N),diag(T)-matrix(1/T,T,T)) %*% matrix(res,N*T,1)) / (sigma_2_eta2*(T-1))
  D22 = ( crossprod(e , tcrossprod(kronecker(Jbar_N,E_T) , et)) ) / (sigma_2_eta2*(T-1))  
  
  D2  = (N/2) * ((1/sigma_2_1)*(D21 - 1) + (T-1) * (1/sigma_2_eta2) * (D22-1) )
  
  LM1 = sqrt(2) * sigma_2_2 * sigma_2_eta1 * (1 / sqrt(T * (T-1) * (sigma_2_eta1^2 + (N-1) * sigma_2_2^2))) * D1 
  LM2 = sqrt(2) * sigma_2_1 * sigma_2_eta2 * (1 / sqrt(N * (N-1) * (sigma_2_eta2^2 + (T-1) * sigma_2_1^2))) * D2
  ### END calculation of parts of test statistic ##

  if (effect=="time") {
    stat = LM1   
    null = c("sigma^2_{gamma}=0 allowing sigma^2_{delta}>=0")
    tested_eff = "individual"
  } else {
    stat = LM2
    null = c("sigma^2_{delta}=0 allowing sigma^2_{gamma}>=0")
    tested_eff = "time"
  }

  pval = pnorm(stat, lower.tail = FALSE)
  
  method = paste("Conditional Lagrange Multiplier Test ", "for ",tested_eff," effects")
  
  RVAL = list(statistic = stat,
                 p.value   = pval,
                 method    = method,
                 null      = null)
  
  class(RVAL) = "lmtest"
  return(RVAL)
}
#data("Grunfeld", package = "plm")
#p = plm(inv ~ value + capital,data = Grunfeld, model = "within")
#clmtest(p)
