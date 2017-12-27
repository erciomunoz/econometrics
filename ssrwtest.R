ssrwtest = function(x,
                    effect = c("individual", "time")) 
{
############## ssrwtest() ############################################
# For a reference see:
# Wooldridge (2003) "Cluster-sample methods in applied econometrics", AER PP 93:133-138
# Written by Ercio Munoz, Graduate Center CUNY
# December 27, 2017
  
  if (x$args$model != "within") {
    x = update(x, model = "within")
  } 
  if (x$args$effect != "twoways") {
    x = update(x, model = "twoways")
  }
 
  pdim = pdim(x)
  N = pdim$nT$n
  T = pdim$nT$T
  balanced = pdim$balanced # true or false

  ### START calculation of parts of test statistic ##
  if (effect=="individual") {
    fx_level = summary(fixef(x, type = "level", effect="individual"))[,c("Estimate","Std. Error")]
    ols      = lm(fx_level[,1]~1,weights=1/fx_level[,2]^2)
    res      = as.matrix(resid(ols))
    SSR      = t(res)%*%res # sum of squared residuals
  } else {
    fx_level = summary(fixef(x, type = "level", effect="time"))[,c("Estimate","Std. Error")]
    ols      = lm(fx_level[,1]~1,weights=1/fx_level[,2]^2)
    res      = as.matrix(resid(ols))
    SSR      = t(res)%*%res # sum of squared residuals 
  }
  ### END calculation of parts of test statistic ##
  
  if (effect=="individual") {
    stat = SSR 
    null = c("sigma^2_{gamma}=0 allowing sigma^2_{delta}>=0")
    tested_eff = "individual"
    G = N
  } else {
    stat = SSR
    null = c("sigma^2_{delta}=0 allowing sigma^2_{gamma}>=0")
    tested_eff = "time"
    G = T
  }
  
  pval = pchisq(stat, df = G-1, lower.tail = FALSE)
  
  method = paste("SSR_w Test ", "for ",tested_eff," effects")
  
  RVAL = list(statistic = stat,
              p.value   = pval,
              method    = method,
              null      = null)
  
  class(RVAL) = "lmtest"
  return(RVAL)
}
