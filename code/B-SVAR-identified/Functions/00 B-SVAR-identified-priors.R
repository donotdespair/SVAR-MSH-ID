

B.SVAR.VW.priors = function(hyperparameters, data, p, M){
   # hyperparameters[1]       = Gamma.lambda.1           # overall shrinkage - own lag
   # hyperparameters[2]       = Gamma.lambda.2           # other variable lag shrinkage
   # hyperparameters[3]       = Gamma.lambda.4           # intercept shrinkage
   # hyperparameters[4]       = Gamma.c                  # lag decay rate
   # hyperparameters[5]       = Gamma.A.VAR              # prior variance of structural matrix
   
   Y              = data
   
   N  = dim(Y)[2]
   T  = dim(Y)[1]
   TT = T - p - 1
   
   Gamma.lambda.1    = hyperparameters[1]
   Gamma.lambda.2    = hyperparameters[2]
   Gamma.lambda.4    = hyperparameters[3]
   Gamma.c           = hyperparameters[4]
   Gamma.A.var       = hyperparameters[5]
   
   A.Var          = vector("list")
   for (i in 1:N){
      A.Var[[i]]  = Gamma.A.var*diag(N)
   }
   
   # Standard deviations of the AR residuals of individual series, 17 lags
   standard.deviations     = array(NA,N)
   for(series in 1:N){
      AR.resid             = ar(Y[,series],aic=FALSE,order.max=17)$resid
      standard.deviations[series]   = sd(AR.resid[!is.na(AR.resid)])
   }
   
   # Litterman prior
   Gamma.sd          = standard.deviations  # Litterman priors 
   
   # Form the variance of the priors
   # Intercepts
   priors.Gamma      = Gamma.sd * Gamma.lambda.4^2
   # Matrix common to all lags
   litterman.sd      = array(1, c(N,N))
   for(j in 1:N){
      litterman.sd[j,]     = Gamma.sd^2 / Gamma.sd[j]^2
   }
   # Loop for lags
   for(lag in 1:p){
      priors.lag           = litterman.sd * Gamma.lambda.1^2 * Gamma.lambda.2^2 * exp( Gamma.c * lag - Gamma.c)
      diag(priors.lag)     = diag(priors.lag) / Gamma.lambda.2^2
      priors.Gamma         = rbind(priors.Gamma, priors.lag)
   }
   
   
   litterman.prior.mean       = array(0,c(1+N*p,N))
   litterman.prior.mean[2:(N+1),] = diag(N)
   litterman.prior.variance   = array(0,c(1+N*p,1+N*p,N))
   for(i in 1:N){
      diag(litterman.prior.variance[,,i]) <- priors.Gamma[,i]
   }
   
   priors   = list(
      Gamma.Var   = litterman.prior.variance,
      Gamma.Mean  = litterman.prior.mean,
      lambda.alpha= 1,
      lambda.beta = 1,
      A.Var       = A.Var,
      w           = matrix(9 * diag(M) + matrix(1,M,M),ncol=1)    # Priors: Hidden Markov Chain, transition probabilities       
   )
   
   return(priors)
}