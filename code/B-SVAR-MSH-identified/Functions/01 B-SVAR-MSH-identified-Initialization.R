#---------------------------------------------------------------------------------------------------
#  B.VECM: Initial value for starting the Gibbs sampler
#---------------------------------------------------------------------------------------------------
  
B.SVAR.MSH.VW.identified.initialization = function(Y, p, M, restrictions, which.Sigma=1){
   # which.Sigma - determines the variance of which variable should be determining the order of lambdas
   data  = as.matrix(Y)
   colnames(data) = colnames(Y)
	N     = dim(Y)[2] 
   T     = dim(Y)[1] 
   TT    = T - p
   d     = sum(restrictions$dj) 
   r  = ncol(restrictions$Q)
	 
   # X, Y 
   #-----------------------------------------------------------------------------
   # Classical VECM estimation
   y.var       = VAR(data, p = p, type = "const")

   U.temp      = residuals(y.var)
   em.output = EM.MSH(Y=U.temp, M=M, p=1, max.iter=1000000, convergence=10e-10, diagonal.PR_TR=0.95, print=TRUE, plot=FALSE) 
   Sigma.temp  = em.output$Sigma
   
   # The first state has the lowest (1,1) element
   ord         = order(Sigma.temp[which.Sigma,which.Sigma,])
   Sigma.temp  = Sigma.temp[,,ord]
   
   lambda.omega      = matrix(NA, N,M)
   lambda.omega[,1]  = as.matrix(diag(Sigma.temp[,,1]))
   for (regime in 2:M){
      lambda.omega[,regime]   = as.vector(diag(Sigma.temp[,,regime])/diag(Sigma.temp[,,1]))
   }
   
   # Structural matrix A (and a)
   a           = as.matrix(rnorm(r,sd=0.01))
   A           = matrix(restrictions$Q%*%a + restrictions$q, ncol=N)
   
   Sigma       = Sigma.temp
   for (regime in 1:M){
      Sigma[,,regime]   = solve(A)%*%diag(as.vector(diag(Sigma.temp[,,regime])))%*%t(solve(A))
   }
   
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

   # PR_TR
   PR_TR       = em.output$PR_TR
   
   hyper       = matrix(1,3,1)
   hyper[2,]   = var(alpha[1,])
   hyper[3,]   = var(matrix(alpha[-1,],ncol=1))
   
   xi          = diag(M)[,c(em.output$Xi[1],em.output$Xi)]
   
	# Output     
   #-----------------------------------------------------------------------------     
   Gibbs.input = list(        
      Y        = data,
      U.SF     = U.SF,
      U        = U,
      Sigma    = Sigma,
      lambda.omega   = lambda.omega,
      A        = A,
      a        = a,
      alpha    = alpha,
      hyper    = hyper,
      PR_TR    = PR_TR,
      w        = matrix(t(PR_TR),ncol=1),
      xi       = xi[ord,],
      MH.logkernel = -Inf
      ) 
   
   return(Gibbs.input)
}
