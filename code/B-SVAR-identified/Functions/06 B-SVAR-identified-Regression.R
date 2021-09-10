
B.SVAR.identified.regression = function(aux, priors){
    
   # Setup constants     
   #-----------------------------------------------------------------------------    
   N  = dim(aux$Y)[2]     
   T  = dim(aux$Y)[1]     
   TT = dim(aux$U)[1]     
   p  = T - TT 
   
   # X, Y 
   #-----------------------------------------------------------------------------
   X        = matrix(data=1, nrow=TT, ncol=1)   # A column of ones for the intercept
   if(p > 0){
      for(lag in 1:p){
         Y.m.j    = aux$Y[(p+1-lag):(T-lag),]
         X        = cbind(X,Y.m.j)
      }
   }
   Y        = matrix(aux$Y[(p+1):T,],ncol=N)   
         
   #-----------------------------------------------------------------------------
   # prior covariance matrix
   prior.cov.inv  = diag(1/c(aux$hyper[2], aux$hyper[3]*diag(priors$beta.H)))
   prior.mean     = rbind(rep(0,N), priors$beta.P)
   
   # lambda.m
   lambda         = aux$lambda.omega

   H     = array(NA, c(1+p*N,1+p*N,N)) 
   P     = vector("list",N)
   for(equation in 1:N){
      XX      = crossprod(X)/lambda[equation,]
      XY      = crossprod(X,Y)/lambda[equation,]
      
      # Parameters of full conditional  
      H[,,equation]     = solve(XX + prior.cov.inv)
      H[,,equation]     = 0.5*(H[,,equation] + t(H[,,equation]))
      P[[equation]]     = H[,,equation] %*% ( XY + prior.cov.inv %*% prior.mean) %*% as.matrix(aux$A[equation,])
   }  
   
   # Draw from multivariate normal (mvtnorm package)
   alpha     = array(NA, c(1+N*p, N))
   for (equation in 1:N) {
      draw              = rmvnorm(n=1, mean=P[[equation]], sigma=H[,,equation], method="chol")
      alpha[,equation]  = draw
   }
        
   # Store    
   #-----------------------------------------------------------------------------    
   # Gamma
   aux$alpha   = alpha
   
   U.SF        = Y %*% t(aux$A) - X %*% alpha
   aux$U.SF    = U.SF
   aux$U       = U.SF %*% solve(t(aux$A))
   
   # Output 
   #-----------------------------------------------------------------------------
   return(aux)
} 
