
B.SVAR.identified.initialization = function(Y, p,restrictions){
   data  = as.matrix(Y)
   colnames(data) = colnames(Y)
	N     = dim(Y)[2] 
   T     = dim(Y)[1] 
   TT    = T - p
   r  = ncol(restrictions$Q)
	 
   # X, Y 
   #-----------------------------------------------------------------------------
   # Classical VECM estimation
   y.var       = VAR(data, p = p, type = "const")

   U.temp      = residuals(y.var)
   Sigma.temp  = crossprod(U.temp)/TT
   
   lambda.omega      = matrix(NA, N,1)
   lambda.omega[,1]  = as.matrix(diag(Sigma.temp))

   # Structural matrix A (and a)
   a           = as.matrix(rnorm(r,sd=0.01))
   A           = matrix(restrictions$Q%*%a + restrictions$q, ncol=N)
   
   # Gamma
   alpha       = t(Bcoef(y.var)[,c(dim(Bcoef(y.var))[2], 1:(dim(Bcoef(y.var))[2]-1))])
   
   # X, Y 
   #-----------------------------------------------------------------------------
   X        = matrix(data=1, nrow=TT, ncol=1)   # A column of ones for the intercept
   if(p > 0){
      for(lag in 1:p){
         Y.m.j    = data[(p+1-lag):(T-lag),]
         X        = cbind(X,Y.m.j)
      }
   }
   Y        = matrix(data[(p+1):T,],ncol=N)   
   
   # U
   U.SF     = Y  %*% t(A) - X %*% alpha
   U        = U.SF %*% solve(t(A))

   hyper       = matrix(1,3,1)
   hyper[2,]   = var(alpha[1,])
   hyper[3,]   = var(matrix(alpha[-1,],ncol=1))
   
   # Output     
   #-----------------------------------------------------------------------------     
   Gibbs.input = list(        
      Y        = data,
      U.SF     = U.SF,
      U        = U,
      Sigma    = Sigma.temp,
      lambda.omega   = lambda.omega,
      A        = A,
      a        = a,
      alpha    = alpha,
      hyper    = hyper,
      MH.logkernel = -Inf
      ) 
   
   return(Gibbs.input)
}
