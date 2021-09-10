
B.SVAR.MSH.identified.structural.MH = function(aux,priors,restrictions,nu=5,C=1){
   N  = dim(aux$Y)[2]     
   T  = dim(aux$Y)[1]     
   TT = dim(aux$U)[1]     
   M  = dim(aux$PR_TR)[1]
   p  = T - TT
   r  = ncol(restrictions$Q)
   
   # X, Y 
   #-----------------------------------------------------------------------------
   X     = matrix(data=1, nrow=TT, ncol=1)   # A column of ones for the intercept
   if(p > 0){
      for(lag in 1:p){
         Y.m.j    = aux$Y[(p+1-lag):(T-lag),]
         X        = cbind(X,Y.m.j)
      }
   }
   Y     = matrix(aux$Y[(p+1):T,],ncol=N)      
   
   T.m      = apply(aux$xi,1,sum)
   lambda.m = aux$lambda.omega
   for (m in 2:M) {
      lambda.m[,m]   = lambda.m[,m]*lambda.m[,1]
   }
   
   X.tilde  = -kronecker(Y,diag(N)) %*% restrictions$Q
   LAMBDA   = diag(1/as.vector(matrix(lambda.m[,apply(aux$xi,2,which.max)],ncol=1)))
   P.tmp    = t(X.tilde) %*% (LAMBDA) %*% X.tilde
   tmp      = try(ginv(P.tmp),silent=FALSE)
#    k        = 1
#    while (inherits(tmp, "try-error") & k<152){
#       P.tmp    = t(X.tilde) %*% (k*LAMBDA) %*% X.tilde
#       tmp      = try(ginv(P.tmp),silent=FALSE)
#       k        = k + 5
#    }
   
   
   if (!inherits(tmp, "try-error")){
      # P.tmp.inv      = k*0.5*(tmp+t(tmp))
      P.tmp.inv      = 0.5*(tmp+t(tmp))
      Y.tilde  = kronecker(Y,diag(N)) %*% restrictions$q -matrix(t(X%*%aux$alpha),ncol=1)
      prodL    = sum(log(lambda.m)%*%diag(-T.m/2))
      
      
      # Current state
      U.tilde.current  = Y.tilde - X.tilde %*% aux$a
      sum.u.current    = as.numeric(t(U.tilde.current) %*% LAMBDA %*% U.tilde.current)
      
      A0.current       = matrix(restrictions$Q %*% aux$a + restrictions$q, ncol=N)
      detA.current     = TT * as.numeric(determinant(A0.current, logarithm=TRUE)$modulus)
      log.likelihood.current = (-TT*N/2)*log(2*pi) + prodL + detA.current - 0.5*sum.u.current
      
      log.prior.a.current    = dmvnorm(as.vector(aux$a), sigma=aux$hyper[1]*diag(r), log=TRUE)
      log.prior.beta.current = 0
      for (n in 1:N){
         log.prior.beta.current = log.prior.beta.current + dmvnorm(aux$alpha[-1,n], mean = A0.current[n,]%*%t(priors$beta.P) , sigma=aux$hyper[3]*diag(p*N), log=TRUE)
      }
      log.kernel.current    = log.likelihood.current + log.prior.a.current + log.prior.beta.current
      
      
      # Candidate
      candidate.a    = t(rmvt(n=1, sigma = C*P.tmp.inv, df = nu, delta = as.vector(aux$a)))
      U.tilde  = Y.tilde - X.tilde %*% candidate.a
      sum.u    = as.numeric(t(U.tilde) %*% LAMBDA %*% U.tilde)
      
      A0       = matrix(restrictions$Q %*% candidate.a + restrictions$q, ncol=N)
      detA     = TT * as.numeric(determinant(A0, logarithm=TRUE)$modulus)
      log.likelihood = (-TT*N/2)*log(2*pi) + prodL + detA - 0.5*sum.u
      
      log.prior.a    = dmvnorm(as.vector(candidate.a), sigma=aux$hyper[1]*diag(r), log=TRUE)
      log.prior.beta = 0
      for (n in 1:N){
         log.prior.beta = log.prior.beta + dmvnorm(aux$alpha[-1,n], mean = A0[n,]%*%t(priors$beta.P) , sigma=aux$hyper[3]*diag(p*N), log=TRUE)
      }
      log.kernel.candidate    = log.likelihood + log.prior.a + log.prior.beta
            
      if (runif(1)< exp(log.kernel.candidate-log.kernel.current)) {
         aux$a    = candidate.a
         aux$A    = A0
         for (regime in 1:M){
            if (regime ==1){
               lambda.m    = aux$lambda.omega[,regime]
            } else {
               lambda.m    = aux$lambda.omega[,regime] * aux$lambda.omega[,1]
            }
            aux$Sigma[,,regime]     = solve(aux$A) %*% diag(lambda.m) %*% solve(t(aux$A))
         }
         aux$MH.logkernel = log.kernel.candidate
      } else {
         aux$MH.logkernel = log.kernel.current
      }
      
   }
   
   return(aux)
}


B.SVAR.MSH.identified.A0.log.kernel = function(a,aux,priors,restrictions){
   N  = dim(aux$Y)[2]     
   T  = dim(aux$Y)[1]     
   TT = dim(aux$U)[1]     
   M  = dim(aux$PR_TR)[1]
   p  = T - TT
   r  = ncol(restrictions$Q)

   # X, Y 
   #-----------------------------------------------------------------------------
   X     = matrix(data=1, nrow=TT, ncol=1)   # A column of ones for the intercept
   if(p > 0){
      for(lag in 1:p){
         Y.m.j    = aux$Y[(p+1-lag):(T-lag),]
         X        = cbind(X,Y.m.j)
      }
   }
   Y     = matrix(aux$Y[(p+1):T,],ncol=N)      
   
   T.m      = apply(aux$xi,1,sum)
   lambda.m = aux$lambda.omega
   for (m in 2:M) {
      lambda.m[,m]   = lambda.m[,m]*lambda.m[,1]
   }
   
   y.tilde  = array(NA,c(N,TT))
   x.tilde  = array(NA,c(N,r,TT))
   u.tilde  = rep(NA,TT)
   for (t in 1:TT){
      y.tilde[,t]    = kronecker(matrix(Y[t,],nrow=1),diag(N)) %*% restrictions$q - t(aux$alpha) %*% X[t,]
      x.tilde[,,t]   = - kronecker(matrix(Y[t,],nrow=1),diag(N)) %*% restrictions$Q
      u.tilde[t]     = t(y.tilde[,t] - x.tilde[,,t] %*% a) %*% diag(1/lambda.m[,which(aux$xi[,t]==1)]) %*% (y.tilde[,t] - x.tilde[,,t] %*% a)
   }
      
   prodL    = sum(log(lambda.m)%*%diag(-T.m/2))
   A0       = matrix(restrictions$Q %*% a + restrictions$q, ncol=N)
   detA     = TT * determinant(A0, logarithm=TRUE)$modulus
   log.likelihood = (-TT*N/2)*log(2*pi) + prodL + detA - 0.5*sum(u.tilde)
   
   log.prior.a    = dmvnorm(as.vector(a), sigma=aux$hyper[1]*diag(r), log=TRUE)
   log.prior.beta = 0
   for (n in 1:N){
      log.prior.beta = log.prior.beta + dmvnorm(aux$alpha[-1,n], mean = A0[n,]%*%t(priors$beta.P) , sigma=aux$hyper[3]*diag(p*N), log=TRUE)
   }
   
   return(as.numeric(log.likelihood + log.prior.a + log.prior.beta))
}
