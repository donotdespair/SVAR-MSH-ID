# Filename  :  04 B-mixture-classification.R
# By        :  tw
# Subject   :  Gibbs step: Sampling classifications for observations

library(abind)
library(fAsianOptions)

ln.d.inverse.gamma.2 = function(x,alpha,beta){
  output = -lgamma(beta/2) + (beta/2)*log(alpha/2) - .5*(beta+2)*log(x) - alpha/(2*x)
  return(output)
}

ln.d.absolute.normal = function(x,mi,ro){
   ln.constant    = ((2*ro+1)/(2*ro))*log(2*ro) + lgamma((1+ro)/(2*ro)) + as.numeric(kummerM(x=-(mu^2)/(2*ro), a=-1/(2*ro), b=.5, lnchf = 1, ip = 0))
   ldan  = (1/ro)*log(abs(x)) - ((x-mi)^2)/(2*ro) - ln.constant
   return(ldan)
}

B.SVAR.MSH.A0.log.kernel = function(A0, theta.star, priors, restrictions){
   # computes the ordinate of a joint kernel of the unrestricted elements of A0
   
   N  = dim(theta.star$Y)[2]     
   T  = dim(theta.star$Y)[1]     
   TT = dim(theta.star$U)[1]     
   p  = T - TT
   
   X     = matrix(data=1, nrow=TT, ncol=1)   # A column of ones for the intercept
   if(p > 0){
      for(lag in 1:p){
         Y.m.j    = theta.star$Y[(p+1-lag):(T-lag),]
         X        = cbind(X,Y.m.j)
      }
   }
   Y     = matrix(theta.star$Y[(p+1):T,],ncol=N)      
   X     = t(X)
   Y     = t(Y)
   
   Omega.bar      = vector("list", N)
   mu.bar         = vector("list", N)
   Omega.bar.tmp  = vector("list", N)
   mu.bar.tmp     = vector("list", N)
   kernel.tmp     = rep(0,1)
   
   for (n in 1:N){
      Omega.bar.tmp[[n]]   = matrix(0, dim(theta.star$Un[[n]])[1], dim(theta.star$Un[[n]])[1])
      mu.bar.tmp[[n]]      = matrix(0, 1, dim(theta.star$Un[[n]])[1])
      for (regime in 1:M){
         xi.m                 = theta.star$xi[regime,]==1
         Omega.bar.tmp[[n]]   = Omega.bar.tmp[[n]] + tcrossprod(theta.star$Un[[n]]%*%Y[,xi.m]/sqrt(theta.star$lambda[n,regime]))
         mu.bar.tmp[[n]]      = mu.bar.tmp[[n]]    + (theta.star$Gamma[,n] %*% X[,xi.m] %*% t(theta.star$Un[[n]] %*% Y[,xi.m])/theta.star$lambda[n,regime])
      }
      Omega.bar[[n]] = solve( Omega.bar.tmp[[n]] + solve(theta.star$Un[[n]] %*% priors$A.Var[[n]] %*% t(theta.star$Un[[n]])) + theta.star$Un[[n]] %*% t(priors$Gamma.Mean) %*% solve(priors$Gamma.Var[,,n]) %*% priors$Gamma.Mean %*% t(theta.star$Un[[n]]))
      mu.bar[[n]]    = (mu.bar.tmp[[n]] + theta.star$Gamma[,n]%*%solve(priors$Gamma.Var[,,n])%*%priors$Gamma.Mean%*%t(theta.star$Un[[n]]) ) %*% Omega.bar[[n]]
      kernel.tmp     = kernel.tmp + (A0[n,diag(restrictions$Q[[n]])==0] - mu.bar[[n]]) %*% solve(Omega.bar[[n]]) %*% t(A0[n,diag(restrictions$Q[[n]])==0] - mu.bar[[n]])
   }

   log.kernel = TT*log(abs(det(A0))) -.5*as.numeric(kernel.tmp)
   return(log.kernel)   
}






B.SVAR.MSH.A0.log.kernel.vectorized = function(A0.array, theta.star, priors, restrictions){
   # The same as B.SVAR.MSH.A0.log.kernel but for a sample of draws of A0 collected in A0.array
   # A vectorized version of B.SVAR.MSH.A0.log.kernel
   
   N  = dim(theta.star$Y)[2]     
   T  = dim(theta.star$Y)[1]     
   TT = dim(theta.star$U)[1]     
   p  = T - TT
   
   X     = matrix(data=1, nrow=TT, ncol=1)   # A column of ones for the intercept
   if(p > 0){
      for(lag in 1:p){
         Y.m.j    = theta.star$Y[(p+1-lag):(T-lag),]
         X        = cbind(X,Y.m.j)
      }
   }
   Y     = matrix(theta.star$Y[(p+1):T,],ncol=N)      
   X     = t(X)
   Y     = t(Y)
   
   Omega.bar      = vector("list", N)
   mu.bar         = vector("list", N)
   Omega.bar.tmp  = vector("list", N)
   mu.bar.tmp     = vector("list", N)
   
   for (n in 1:N){
      Omega.bar.tmp[[n]]   = matrix(0, dim(theta.star$Un[[n]])[1], dim(theta.star$Un[[n]])[1])
      mu.bar.tmp[[n]]      = matrix(0, 1, dim(theta.star$Un[[n]])[1])
      for (regime in 1:M){
         xi.m                 = theta.star$xi[regime,]==1
         Omega.bar.tmp[[n]]   = Omega.bar.tmp[[n]] + tcrossprod(theta.star$Un[[n]]%*%Y[,xi.m]/sqrt(theta.star$lambda[n,regime]))
         mu.bar.tmp[[n]]      = mu.bar.tmp[[n]]    + (theta.star$Gamma[,n] %*% X[,xi.m] %*% t(theta.star$Un[[n]] %*% Y[,xi.m])/theta.star$lambda[n,regime])
      }
      Omega.bar[[n]] = solve( Omega.bar.tmp[[n]] + solve(theta.star$Un[[n]] %*% priors$A.Var[[n]] %*% t(theta.star$Un[[n]])) + theta.star$Un[[n]] %*% t(priors$Gamma.Mean) %*% solve(priors$Gamma.Var[,,n]) %*% priors$Gamma.Mean %*% t(theta.star$Un[[n]]))
      mu.bar[[n]]    = (mu.bar.tmp[[n]] + theta.star$Gamma[,n]%*%solve(priors$Gamma.Var[,,n])%*%priors$Gamma.Mean%*%t(theta.star$Un[[n]]) ) %*% Omega.bar[[n]]
   }
   
   log.kernel.A0 = function(A0){
      kernel.tmp     = rep(0,1)
      for (n in 1:N){
         kernel.tmp     = kernel.tmp + (A0[n,diag(restrictions$Q[[n]])==0] - mu.bar[[n]]) %*% solve(Omega.bar[[n]]) %*% t(A0[n,diag(restrictions$Q[[n]])==0] - mu.bar[[n]])
      }
      log.kernel = TT*log(abs(det(A0))) -.5*as.numeric(kernel.tmp)
      return(log.kernel)  
   }
   
   log.kernel = apply(A0.array,3,log.kernel.A0)
   return(log.kernel)   
}






B.SVAR.MSH.VW.draw.full.conditional.A0 = function(S, theta.star, priors, restrictions){
   N  = dim(theta.star$Y)[2]     
   T  = dim(theta.star$Y)[1]     
   TT = dim(theta.star$U)[1]     
   p  = T - TT
   
   X     = matrix(data=1, nrow=TT, ncol=1)   # A column of ones for the intercept
   if(p > 0){
      for(lag in 1:p){
         Y.m.j    = theta.star$Y[(p+1-lag):(T-lag),]
         X        = cbind(X,Y.m.j)
      }
   }
   Y     = matrix(theta.star$Y[(p+1):T,],ncol=N)      
   X     = t(X)
   Y     = t(Y)
   
   Omega.bar      = vector("list", N)
   mu.bar         = vector("list", N)
   Omega.bar.tmp  = vector("list", N)
   mu.bar.tmp     = vector("list", N)
   
   for (n in 1:N){
      Omega.bar.tmp[[n]]   = matrix(0, dim(theta.star$Un[[n]])[1], dim(theta.star$Un[[n]])[1])
      mu.bar.tmp[[n]]      = matrix(0, 1, dim(theta.star$Un[[n]])[1])
      for (regime in 1:M){
         xi.m                 = theta.star$xi[regime,]==1
         Omega.bar.tmp[[n]]   = Omega.bar.tmp[[n]] + tcrossprod(theta.star$Un[[n]]%*%Y[,xi.m]/sqrt(theta.star$lambda[n,regime]))
         mu.bar.tmp[[n]]      = mu.bar.tmp[[n]]    + (theta.star$Gamma[,n] %*% X[,xi.m] %*% t(theta.star$Un[[n]] %*% Y[,xi.m])/theta.star$lambda[n,regime])
      }
      Omega.bar[[n]] = solve( Omega.bar.tmp[[n]] + solve(theta.star$Un[[n]] %*% priors$A.Var[[n]] %*% t(theta.star$Un[[n]])) + theta.star$Un[[n]] %*% t(priors$Gamma.Mean) %*% solve(priors$Gamma.Var[,,n]) %*% priors$Gamma.Mean %*% t(theta.star$Un[[n]]))
      mu.bar[[n]]    = (mu.bar.tmp[[n]] + theta.star$Gamma[,n]%*%solve(priors$Gamma.Var[,,n])%*%priors$Gamma.Mean%*%t(theta.star$Un[[n]]) ) %*% Omega.bar[[n]]
   }
   
   A0.tmp         = theta.star$A
   A0.array       = array(NA,c(dim(theta.star$A),S))
   for (s in 1:S){
      w     = array(0, c(N,1))
      for(i.star in 1:N){     # given last A0gbs and generate new A0bgs  
         A.star   = t(A0.tmp[-i.star,])
         w        = orthogonal.complement.matrix.TW(A.star)
         
         #*** Constructing orthonormal basis w_1, ..., w_qi at each Gibbs step
         Tn    = t(chol(TT*Omega.bar[[i.star]]))    # defined below eq (14)
         w0    = t(Tn) %*% theta.star$Un[[i.star]] %*% w
         w1    = w0 / as.numeric(sqrt(crossprod(w0)))
         if (dim(theta.star$Un[[i.star]])[1]>1) {
            W     = orthogonal.complement.matrix.TW(w1)
         }
         
         # compute the means of the Absolute Normal and normal distributions distribution:
         xi    = rep(NA,dim(theta.star$Un[[i.star]])[1])
         xi[1] = mu.bar[[i.star]] %*% solve(t(Tn)) %*% w1
         if (length(xi)>1){
            for (j in 2:length(xi)){
               xi[j]    = mu.bar[[i.star]] %*% solve(t(Tn)) %*% W[,j-1]
            }
         }
         
         #*** Draw beta's according to Proposition C.1 of Villani(2009, JAE).
         # All but first element
         gkbeta         = array(NA, c(dim(theta.star$Un[[i.star]])[1],1))  # qi-by-1: greak beta's
         if (length(xi)>1){
            for (j in 2:length(xi)){
               gkbeta[j,]    = rnorm(n=1, mean=xi[j], sd=sqrt(1/TT))
            }
         }
         
         # Sample first element beta1
         # A draw from the absolute normal distribution according to 
         mu.an       = rep(NA,2)
         s.an        = rep(NA,2)
         mu.an[1]    = .5*xi[1] - .5*sqrt(xi[1]^2+4)
         mu.an[2]    = .5*xi[1] + .5*sqrt(xi[1]^2+4)
         for (j in 1:2) s.an[j]  = ((mu.an[j]^2)*(1/TT))/(1+mu.an[j]^2)
         www         = 1/(1+exp(2*xi[1]*TT))
         param       = norMix(mu=mu.an, sigma = sqrt(s.an), w = c(www,1-www))
         gkbeta[1,]  = rnorMix(n=1, obj=param)
         
         if (length(xi)>1) {
            WWW = cbind(w1,W)
         } else {
            WWW = w1
         }
         A0.array[i.star,,s]    = as.vector(Tn %*% WWW %*% gkbeta) %*% theta.star$Un[[i.star]]
      }
   }
   
   return(A0.array)
}



B.SVAR.MSH.A0.log.q.vectorized = function(S, A0.array){
   A0.mean     = apply(A0.array,1:2,mean)
   
   A0.q.draw   = array(0,c(dim(A0.array)[1:2],S))
   dg.tmp      = matrix(NA,dim(A0.array)[1],S)
   q.mean      = vector("list",dim(A0.array)[1])
   q.cov       = vector("list",dim(A0.array)[1])
   for (n in 1:dim(A0.array)[1]){
      non.zero    = A0.mean[n,]!=0
      A0n.mean    = A0.mean[n,non.zero]
      if (sum(non.zero)==1){
         A0n.cov     = cov(as.matrix(A0.array[n,non.zero,]))
         A0.q.draw[n,non.zero,] = t(rmvnorm(n=S, mean = A0n.mean, sigma = A0n.cov))
         dg.tmp[n,]  = apply(matrix(A0.q.draw[n,non.zero,],nrow=1),2, dmvnorm, mean = A0n.mean, sigma = A0n.cov, log=TRUE)
      } else {
         A0n.cov     = cov(t(as.matrix(A0.array[n,non.zero,])))
         A0.q.draw[n,non.zero,] = t(rmvnorm(n=S, mean = A0n.mean, sigma = A0n.cov))
         dg.tmp[n,]  = apply(A0.q.draw[n,non.zero,],2, dmvnorm, mean = A0n.mean, sigma = A0n.cov, log=TRUE)
      }
      q.mean[[n]] = A0n.mean
      q.cov[[n]] = A0n.cov      
   }
   
   output         = new.env()
   output$draws   = A0.q.draw
   output$log.pdf = apply(dg.tmp,2,sum)
   output$q.mean  = q.mean
   output$q.cov   = q.cov
   return(as.list(output))
}



B.SVAR.MSH.A0.kernel.normalizing.constant = function(S, theta.star, priors, restrictions){
   A0.full.draws     = B.SVAR.MSH.VW.draw.full.conditional.A0(S, theta.star, priors, restrictions)
   A0.log.kernel     = B.SVAR.MSH.A0.log.kernel.vectorized(A0.array=A0.full.draws, theta.star, priors, restrictions)
   A0.q.draws        = B.SVAR.MSH.A0.log.q.vectorized(S, A0.array=A0.full.draws)
   A0.q.log.kernel   = B.SVAR.MSH.A0.log.kernel.vectorized(A0.array=A0.q.draws$draws, theta.star, priors, restrictions)
   logc              = mean(A0.q.draws$log.pdf - A0.log.kernel)
   log.mdd           = logc + log(sum(exp(A0.q.draws$log.pdf - A0.log.kernel - logc))) - log(S)
   return(log.mdd)
}


B.SVAR.MSH.A0.log.pdf.vectorized = function(A0.array, theta.star, priors, restrictions, S.constant=1000){
   
   if (is.matrix(A0.array)) {
      A0.tmp         = array(NA,c(dim(A0.array),1))
      A0.tmp[,,1]    = A0.array
      A0.array       = A0.tmp
   }
   
   log.k       = B.SVAR.MSH.A0.log.kernel.vectorized(A0.array, theta.star, priors, restrictions)
   log.c       = B.SVAR.MSH.A0.kernel.normalizing.constant(S=S.constant, theta.star, priors, restrictions)
   log.d       = log.k + log.c
   return(log.d)   
}







MC.hist = function(st){
   st.S           = dim(st)[2]
   st.unique      = matrix(NA,dim(st)[1],0)
   st.distance    = function(S1, S2){ length(S1) - sum(S1==S2) }
   st.which       = vector("list")
   st.count       = c()
   indicator      = 1
   while (ncol(st)!=0){
      st.dist     = apply(st,2,st.distance, S2=st[,1])
      st.which[[indicator]]   = which(st.dist==0)
      st.unique   = cbind(st.unique,st[,1])
      st.count    = c(st.count,length(st.which[[indicator]]))
      st          = as.matrix(st[,-st.which[[indicator]]])
      indicator   = indicator + 1
   }
   st.freq        = st.count/st.S
   st.order       = order(st.count, decreasing=TRUE)
   st.unique      = st.unique[,st.order]
   st.freq        = st.freq[st.order]
   st.count       = st.count[st.order]
   st.which.tmp   = vector("list")
   for (i in 1:length(st.order)){
      st.which.tmp[[i]]    = st.which[[st.order[i]]]
   }
   output         = new.env()
   output$unique  = as.matrix(st.unique)
   output$freq    = st.freq
   output$count   = st.count
   output$which   = st.which.tmp
   return(as.list(output))
}








B.SVAR.MSH.S.bar = function(S.t){
   S     = dim(S.t)[2]
   M     = max(S.t)
   TT    = dim(S.t)[1]
   
   xi    = array(NA,c(M,TT,S))
   for (s in 1:S){
      xi[,,s]    = diag(M)[,S.t[,s]]
   }
   posterior= apply(xi,1:2,mean)
#    plot.ts(t(posterior))
   S.bar    = apply(posterior,2,which.max)
   return(S.bar)
}

#######################################################
# PR_TR full conditional posterior ordinate VECTORIZED
#######################################################

B.SVAR.MSH.PR_TR.log.pdf = function(S.t, theta.star.g, priors, restrictions){
   d     = sum(restrictions$dj) 
   Q     = length(restrictions$dj)
   
   transitions.matrix   = count.regime.transitions(diag(M)[,S.t])
   transitions.vector   = matrix(transitions.matrix, ncol=1)
   wj.transitions       = array(NA, c(d,1)) 
   log.full             = rep(NA,Q)
   for(j in 1:Q){
      dj                = restrictions$dj[j]
      if (dj>1){
         index.begin    = sum(restrictions$dj[1:j-1]) + 1
         index.end      = index.begin + dj - 1
         for(i in index.begin:index.end){
            wj.transitions[i]       = matrix(restrictions$M[,i],nrow=1) %*% transitions.vector 
         }
         alpha          = wj.transitions[index.begin:index.end] + priors$w[index.begin:index.end]
         log.full[j]    = log(ddirichlet(x=theta.star.g$w[index.begin:index.end], alpha=alpha))
      }
   }
   return(sum(log.full))
}
   











B.SVAR.MSH.mdd = function(theta.star, qqq, priors, restrictions, permutation.fraction, S.aux, S.constant, keep.every=1, min.count.no=5, debug=FALSE){
   
   cat("\nEstimation of the MDD for the SVAR-MSH model")
   
   S  = dim(qqq$posteriors$lambda)[3]
   N  = dim(qqq$last.draws$Y)[2]
   M  = dim(qqq$last.draws$PR_TR)[1]
   p  = (dim(theta.star$Y)[1] - dim(theta.star$U)[1])
   T  = dim(theta.star$Y)[1]
   TT = T - p
   d     = sum(restrictions$dj) 
   Q     = length(restrictions$dj)
   
   rest  = 0
   for (n in 1:N){
      rest           = rest + sum(diag(restrictions$Q[[n]])==1)
   }
   if (rest!=0){
      restricted     = TRUE
      permutation.fraction    = 1
      cat(" with restricted matrix of contemporaneous effects, A0")
   } else {
      restricted     = FALSE
      cat("\nwith unrestricted matrix of contemporaneous effects, A0")
   }
   
   # log of the lieklihood function ordinate at theta.star
   eta.t             = matrix(NA, M, TT)
   for(regime in 1:M){
      eta.tmp        = dmvnorm(theta.star$U, sigma=theta.star$Sigma[,,regime]) 
      eta.t[regime,] = matrix(eta.tmp, nrow=1)
   }
   PR.forecasted     = G.BHLK.forecasting(theta.star)
   log.likelihood    = sum( log(colSums(eta.t * PR.forecasted)))
   
   # log of the prior pdf ordinate at theta.star
   log.prior.Gamma   = rep(NA,N)
   for (n in 1:N){
      log.prior.Gamma[n]   = dmvnorm(x=theta.star$Gamma[,n] , mean=priors$Gamma.Mean %*% theta.star$A[n,], sigma=priors$Gamma.Var[,,n], log=TRUE)
   }
   log.prior.A       = rep(NA,N)
   for (n in 1:N){
      log.prior.A[n]       = dmvnorm(x=theta.star$a[[n]] , mean=rep(0,length(theta.star$a[[n]])), sigma=theta.star$Un[[n]] %*% priors$A.Var[[n]] %*% t(theta.star$Un[[n]]), log=TRUE)
   }
   log.prior.lambda   = matrix(NA,N,M-1)
   for (n in 1:N){
      for (m in 2:M){
         log.prior.lambda[n,m-1]    = ln.d.inverse.gamma.2(x=theta.star$lambda[n,m], alpha=priors$lambda.alpha, beta=priors$lambda.beta)
      }
   }
#    log.prior.PR_TR   = rep(NA,M)
#    for (m in 1:M){
#       log.prior.PR_TR[m]            = log(ddirichlet(x=theta.star$PR_TR[m,], alpha=priors$w[m,]))
#    }
   
   log.prior.PR_TR   = rep(NA,Q)
   for(j in 1:Q){
      dj                = restrictions$dj[j]
      if (dj>1){
         index.begin    = sum(restrictions$dj[1:j-1]) + 1
         index.end      = index.begin + dj - 1
         log.prior.PR_TR[j]    = log(ddirichlet(x=theta.star$w[index.begin:index.end], alpha=priors$w[index.begin:index.end]))
      }
   }
   
   log.prior         = sum(log.prior.Gamma) + sum(log.prior.A) + sum(log.prior.lambda) + sum(log.prior.PR_TR)
   
   cat("\nValue of the log-likelihood function: ",log.likelihood)
   cat("\nValue of the log of the prior distribution: ",log.prior)

   # log of the posterior density ordinate at theta.star
   # step 1: compute permutations
   ##################################################################
   theta.star.p      = B.SVAR.MSH.UNRESTRICTED.state.labels.permutation.augmentation(theta.star)
   if (!restricted){
      theta.star.p      = B.SVAR.MSH.UNRESTRICTED.equation.ordering.permutation.augmentation(theta.star.p)
   }
   permutations.no   = floor(dim(theta.star.p$lambda)[3]*permutation.fraction)
   if (!restricted){
      perm.indexes      = sample.int(n=dim(theta.star.p$lambda)[3], size = permutations.no)
      cat("\nEstimator computed on the basis of ",length(perm.indexes), "random permutations:\n")
   } else {
      perm.indexes      = 1:dim(theta.star.p$lambda)[3]
      cat("\nEstimator computed on the basis of a complete set of ",length(perm.indexes), " permutations:\n")
   }
   
   # step 2: ordinate of alpha.star | y
   ##################################################################
   # X, Y 
   X        = matrix(data=1, nrow=TT, ncol=1)   # A column of ones for the intercept
   if(p > 0){
      for(lag in 1:p){
         Y.m.j    = theta.star$Y[(p+1-lag):(T-lag),]
         X        = cbind(X,Y.m.j)
      }
   }
   Y        = matrix(theta.star$Y[(p+1):T,],ncol=N)   


   # Step 3.1: log-ordinate of alpha.star | y
   # this part uses the posterior draws
   log.full.Gamma.s     = matrix(NA,N,S)
   H                    = array(NA, dim(priors$Gamma.Var)) 
   P                    = vector("list",N)
   for (s in 1:S){   
      for(equation in 1:N){
         # Crossproducts
         XX          = matrix(0, ncol(X), ncol(X))
         XY          = matrix(0, ncol(X), N)      
         for (regime in 1:M){
            xi.m    = qqq$posteriors$S.t[,s] == regime
            XX      = XX + crossprod(X[xi.m,])/theta.star$lambda[equation,regime]
            XY      = XY + crossprod(X[xi.m,],Y[xi.m,])/theta.star$lambda[equation,regime]
         }
         H[,,equation]  = solve(XX + solve(priors$Gamma.Var[,,equation]))
         P[[equation]]  = H[,,equation] %*% ( XY + solve(priors$Gamma.Var[,,equation]) %*% priors$Gamma.Mean) %*% as.matrix(theta.star$A[equation,])
         log.full.Gamma.s[equation,s]     = dmvnorm(x=theta.star$Gamma[,equation], mean=P[[equation]], sigma=H[,,equation], log=TRUE)
      }  
   }
   
   theta.star.g         = theta.star
#    log.full.Gamma.s     = array(NA,c(N,S,length(perm.indexes)))
   log.full.lambda.s    = array(NA,c((S.aux/keep.every),length(perm.indexes) ))
   log.full.A0.s        = c()
   log.full.PR_TR.s     = array(NA,c((S.aux/keep.every),length(perm.indexes) ))

   t0    = proc.time()
   for (g in 1:length(perm.indexes)){
      if (debug & g%%1==0) {cat("\nPermutation No.: ",g)}
      theta.star.g$Sigma      = theta.star.p$Sigma[,,,perm.indexes[g]]
      theta.star.g$lambda     = theta.star.p$lambda[,,perm.indexes[g]]
      theta.star.g$A          = theta.star.p$A[,,perm.indexes[g]]
      for (n in 1:N) {theta.star.g$a[[n]]     = theta.star.p$a[[n]][,perm.indexes[g]]}
      theta.star.g$Gamma      = theta.star.p$Gamma[,,perm.indexes[g]]
      theta.star.g$PR_TR      = theta.star.p$PR_TR[,,perm.indexes[g]]
      theta.star.g$w          = theta.star.p$w[,perm.indexes[g]]
      theta.star.g$xi         = theta.star.p$xi[,,perm.indexes[g]]
      
#       # Step 3.1: log-ordinate of alpha.star | y
#       # this part uses the posterior draws
#       H     = array(NA, dim(priors$Gamma.Var)) 
#       P     = vector("list",N)
#       label.perm        = allPermutations(M)[perm.indexes[g]%/%prod(1:N) + 1,]
#       row.perm          = allPermutations(N)[perm.indexes[g]%%prod(1:N) ,]
#       for (s in 1:S){
#          # Recreate A0 and lambda after label and row permutations
#          A0       = qqq$posteriors$A[,,s]
#          la       = qqq$posteriors$lambda[,,s]
#          sign.A0  = sign(diag(A0))
#          A0.tmp   = diag(sign.A0)%*%A0
#          tau1     = 1/diag(A0.tmp)^2
#          B0       = diag(sqrt(tau1))%*%A0.tmp
#          tau = la
#          for (m in 1:M) {
#             tau[,m] = la[,m]*tau1
#          }
#          tau.l    = tau[,label.perm]
#          A0.l     = diag(1/sqrt(tau.l[,1]))%*%B0
#          A0.l     = diag(sign.A0)%*%A0.l
#          lambda.l = tau.l
#          for (m in 1:M) {
#             lambda.l[,m] = tau.l[,m]/tau.l[,1]
#          }
#          A0.l     = A0.l[row.perm,]
#          lambda.l = lambda.l[row.perm,]
#          
#          for(equation in 1:N){
#             # Crossproducts
#             XX   = matrix(0, ncol(X), ncol(X))
#             XY   = matrix(0, ncol(X), N)      
#             for (regime in 1:M){
#                xi.m    = qqq$posteriors$S.t[,s] == regime
#                XX      = XX + crossprod(X[xi.m,])/lambda.l[equation,regime]
#                XY      = XY + crossprod(X[xi.m,],Y[xi.m,])/lambda.l[equation,regime]
#             }
#             H[,,equation]     = solve(XX + solve(priors$Gamma.Var[,,equation]))
#             P[[equation]]     = H[,,equation] %*% ( XY + solve(priors$Gamma.Var[,,equation]) %*% priors$Gamma.Mean) %*% as.matrix(A0.l[equation,])
#             log.full.Gamma.s[equation,s,g] = dmvnorm(x=theta.star.p$Gamma[,equation,g], mean=P[[equation]], sigma=H[,,equation], log=TRUE)
#          }  
#       }
      
      # step 3.2: ordinate of lambda.star | y, alpha.star 
      qqq.g    = B.SVAR.MSH.reduced.sampler.log.full.lambda(S.aux, priors, restrictions, starting.values=theta.star.g, debug=FALSE)
      qqq.g    = B.SVAR.MSH.short.mcmc(qqq.g,keep.every)
      for (ss in 1:(S.aux/keep.every)){   
         U.SF           = Y %*% t(qqq.g$posteriors$A[,,ss]) - X %*% theta.star.g$Gamma
         lambda.beta.full.tmp    = matrix(NA,N,M-1)
         for(regime in 2:M){
            # Crossproduct of the residuals
            xi.m              = qqq.g$posteriors$S.t[,ss] == regime
            T.m               = sum(xi.m)
            lambda.beta.full  = priors$lambda.beta + T.m
            if (T.m != 0) {
               Sigma.epsilon  = crossprod(U.SF[xi.m,])
               for(equation in 1:N){
                  lambda.alpha.full       = priors$lambda.alpha + Sigma.epsilon[equation,equation]
                  lambda.beta.full.tmp[equation,regime-1]    = ln.d.inverse.gamma.2(x=theta.star.g$lambda[equation,regime], alpha=lambda.alpha.full, beta=lambda.beta.full)
               }
            }
         }
         log.full.lambda.s[ss,g]       = sum(lambda.beta.full.tmp)
      }
      
      # step 3.3: ordinate of A.star | y, alpha.star, lambda.star
      qqq.g    = B.SVAR.MSH.reduced.sampler.log.full.an(S.aux, priors, restrictions, starting.values=theta.star.g, debug=FALSE)
      qqq.g    = B.SVAR.MSH.short.mcmc(qqq.g,keep.every)
      st.hist  = MC.hist(qqq.g$posteriors$S.t)
      theta.star.tmp       = theta.star.g
      for (ind in 1:max(which(st.hist$count>min.count.no))){
         theta.star.tmp$xi = diag(M)[,st.hist$unique[,ind]]
         log.full.A0.s     = c(as.vector(log.full.A0.s),B.SVAR.MSH.A0.log.pdf.vectorized(A0.array = qqq.g$posteriors$A[,,st.hist$which[[ind]]], theta.star = theta.star.tmp, priors, restrictions, S.constant=S.constant))
      }
      
      # step 3.4: ordinate of P.star | y, alpha.star, lambda.star, A.star
      qqq.g    = B.SVAR.MSH.reduced.sampler.log.full.PR_TR(S.aux, priors, restrictions, starting.values=theta.star.g, debug=FALSE)
      qqq.g    = B.SVAR.MSH.short.mcmc(qqq.g,keep.every)
      log.full.PR_TR.s[,g]        = apply(qqq.g$posteriors$S.t,2,B.SVAR.MSH.PR_TR.log.pdf, theta.star.g=theta.star.g, priors=priors, restrictions=restrictions)   
      t1    = proc.time()
      if (debug & g%%1==0) {cat(", total time elapsed: ", (t1-t0)[1]/60," [min]")}
   }

   lnc.lambda           = quantile(log.full.lambda.s, probs=seq(from=.5, to=.9999999, by=.001))
   log.full.lambda.p    = rep(NA,length(lnc.lambda))
   for (i in 1:length(lnc.lambda)){
      log.full.lambda.p[i]    = lnc.lambda[i] + log(mean(exp(log.full.lambda.s - lnc.lambda[i])))
   }
   log.full.lambda      = log.full.lambda.p[which(is.finite(log.full.lambda.p))[1]]

   log.full.Gamma.s.tmp = apply(log.full.Gamma.s,2,sum)
   lnc.Gamma            = quantile(log.full.Gamma.s, probs=seq(from=.5, to=.9999999, by=.001))
   log.full.Gamma.p     = rep(NA,length(lnc.Gamma))
   for (i in 1:length(lnc.Gamma)){
      log.full.Gamma.p[i]  = lnc.Gamma[i] + log(mean(exp(log.full.Gamma.s.tmp - lnc.Gamma[i])))
   }
   log.full.Gamma       = log.full.Gamma.p[which(is.finite(log.full.Gamma.p))[1]]

   lnc.A0               = quantile(log.full.A0.s, probs=seq(from=.5, to=.9999999, by=.0001))
   log.full.A0.p        = rep(NA,length(lnc.A0))
   for (i in 1:length(lnc.A0)){
      log.full.A0.p[i]  = lnc.A0[i] + log(mean(exp(log.full.A0.s - lnc.A0[i])))
   }
   log.full.A0          = log.full.A0.p[which(is.finite(log.full.A0.p))[1]]

   log.full.PR_TR       = mean(na.omit(log.full.PR_TR.s)) + log(mean(exp(na.omit(log.full.PR_TR.s) - mean(na.omit(log.full.PR_TR.s)))))

   # final computation
   lnmdd = log.likelihood + log.prior - log.full.Gamma - log.full.lambda - log.full.A0 - log.full.PR_TR
   cat("\nlog-MDD value: ",lnmdd)
   output         = new.env()
   output$lnmdd   = lnmdd
   output$components    = list(
         log.likelihood = log.likelihood,
         log.prior      = log.prior,
         log.full.Gamma = log.full.Gamma,
         log.full.lambda= log.full.lambda,
         log.full.A0    = log.full.A0,
         log.full.PR_TR = log.full.PR_TR
      )
   output$lnpdf         = list(
         log.full.Gamma.s  = log.full.Gamma.s.tmp,
         log.full.lambda.s = log.full.lambda.s,
         log.full.A0.s     = log.full.A0.s,
         log.full.PR_TR.s  = log.full.PR_TR.s
      )

   return(as.list(output))
}




























B.SVAR.MSH.mdd.justLambda = function(theta.star, qqq, priors, restrictions, permutation.fraction, S.aux, S.constant, keep.every=1, min.count.no=5, debug=FALSE){
  
  cat("\nEstimation of the MDD for the SVAR-MSH model")
  
  S  = dim(qqq$posteriors$lambda)[3]
  N  = dim(qqq$last.draws$Y)[2]
  M  = dim(qqq$last.draws$PR_TR)[1]
  p  = (dim(theta.star$Y)[1] - dim(theta.star$U)[1])
  T  = dim(theta.star$Y)[1]
  TT = T - p
  d     = sum(restrictions$dj) 
  Q     = length(restrictions$dj)
  
  rest  = 0
  for (n in 1:N){
    rest           = rest + sum(diag(restrictions$Q[[n]])==1)
  }
  if (rest!=0){
    restricted     = TRUE
    permutation.fraction    = 1
    cat(" with restricted matrix of contemporaneous effects, A0")
  } else {
    restricted     = FALSE
    cat("\nwith unrestricted matrix of contemporaneous effects, A0")
  }
  
#   # log of the lieklihood function ordinate at theta.star
#   eta.t             = matrix(NA, M, TT)
#   for(regime in 1:M){
#     eta.tmp        = dmvnorm(theta.star$U, sigma=theta.star$Sigma[,,regime]) 
#     eta.t[regime,] = matrix(eta.tmp, nrow=1)
#   }
#   PR.forecasted     = G.BHLK.forecasting(theta.star)
#   log.likelihood    = sum( log(colSums(eta.t * PR.forecasted)))
#   
#   # log of the prior pdf ordinate at theta.star
#   log.prior.Gamma   = rep(NA,N)
#   for (n in 1:N){
#     log.prior.Gamma[n]   = dmvnorm(x=theta.star$Gamma[,n] , mean=priors$Gamma.Mean %*% theta.star$A[n,], sigma=priors$Gamma.Var[,,n], log=TRUE)
#   }
#   log.prior.A       = rep(NA,N)
#   for (n in 1:N){
#     log.prior.A[n]       = dmvnorm(x=theta.star$a[[n]] , mean=rep(0,length(theta.star$a[[n]])), sigma=theta.star$Un[[n]] %*% priors$A.Var[[n]] %*% t(theta.star$Un[[n]]), log=TRUE)
#   }
#   log.prior.lambda   = matrix(NA,N,M-1)
#   for (n in 1:N){
#     for (m in 2:M){
#       log.prior.lambda[n,m-1]    = ln.d.inverse.gamma.2(x=theta.star$lambda[n,m], alpha=priors$lambda.alpha, beta=priors$lambda.beta)
#     }
#   }
#   #    log.prior.PR_TR   = rep(NA,M)
#   #    for (m in 1:M){
#   #       log.prior.PR_TR[m]            = log(ddirichlet(x=theta.star$PR_TR[m,], alpha=priors$w[m,]))
#   #    }
#   
#   log.prior.PR_TR   = rep(NA,Q)
#   for(j in 1:Q){
#     dj                = restrictions$dj[j]
#     if (dj>1){
#       index.begin    = sum(restrictions$dj[1:j-1]) + 1
#       index.end      = index.begin + dj - 1
#       log.prior.PR_TR[j]    = log(ddirichlet(x=theta.star$w[index.begin:index.end], alpha=priors$w[index.begin:index.end]))
#     }
#   }
#   
#   log.prior         = sum(log.prior.Gamma) + sum(log.prior.A) + sum(log.prior.lambda) + sum(log.prior.PR_TR)
#   
#   cat("\nValue of the log-likelihood function: ",log.likelihood)
#   cat("\nValue of the log of the prior distribution: ",log.prior)
  
  # log of the posterior density ordinate at theta.star
  # step 1: compute permutations
  ##################################################################
  theta.star.p      = B.SVAR.MSH.UNRESTRICTED.state.labels.permutation.augmentation(theta.star)
  if (!restricted){
    theta.star.p      = B.SVAR.MSH.UNRESTRICTED.equation.ordering.permutation.augmentation(theta.star.p)
  }
  permutations.no   = floor(dim(theta.star.p$lambda)[3]*permutation.fraction)
  if (!restricted){
    perm.indexes      = sample.int(n=dim(theta.star.p$lambda)[3], size = permutations.no)
    cat("\nEstimator computed on the basis of ",length(perm.indexes), "random permutations:\n")
  } else {
    perm.indexes      = 1:dim(theta.star.p$lambda)[3]
    cat("\nEstimator computed on the basis of a complete set of ",length(perm.indexes), " permutations:\n")
  }
  
  # step 2: ordinate of alpha.star | y
  ##################################################################
  # X, Y 
  X        = matrix(data=1, nrow=TT, ncol=1)   # A column of ones for the intercept
  if(p > 0){
    for(lag in 1:p){
      Y.m.j    = theta.star$Y[(p+1-lag):(T-lag),]
      X        = cbind(X,Y.m.j)
    }
  }
  Y        = matrix(theta.star$Y[(p+1):T,],ncol=N)   
  
  
#   # Step 3.1: log-ordinate of alpha.star | y
#   # this part uses the posterior draws
#   log.full.Gamma.s     = matrix(NA,N,S)
#   H                    = array(NA, dim(priors$Gamma.Var)) 
#   P                    = vector("list",N)
#   for (s in 1:S){   
#     for(equation in 1:N){
#       # Crossproducts
#       XX          = matrix(0, ncol(X), ncol(X))
#       XY          = matrix(0, ncol(X), N)      
#       for (regime in 1:M){
#         xi.m    = qqq$posteriors$S.t[,s] == regime
#         XX      = XX + crossprod(X[xi.m,])/theta.star$lambda[equation,regime]
#         XY      = XY + crossprod(X[xi.m,],Y[xi.m,])/theta.star$lambda[equation,regime]
#       }
#       H[,,equation]  = solve(XX + solve(priors$Gamma.Var[,,equation]))
#       P[[equation]]  = H[,,equation] %*% ( XY + solve(priors$Gamma.Var[,,equation]) %*% priors$Gamma.Mean) %*% as.matrix(theta.star$A[equation,])
#       log.full.Gamma.s[equation,s]     = dmvnorm(x=theta.star$Gamma[,equation], mean=P[[equation]], sigma=H[,,equation], log=TRUE)
#     }  
#   }
  
  theta.star.g         = theta.star
  #    log.full.Gamma.s     = array(NA,c(N,S,length(perm.indexes)))
  log.full.lambda.s    = array(NA,c((S.aux/keep.every),length(perm.indexes) ))
#   log.full.A0.s        = c()
#   log.full.PR_TR.s     = array(NA,c((S.aux/keep.every),length(perm.indexes) ))
  
  t0    = proc.time()
  for (g in 1:length(perm.indexes)){
    if (debug & g%%1==0) {cat("\nPermutation No.: ",g)}
    theta.star.g$Sigma      = theta.star.p$Sigma[,,,perm.indexes[g]]
    theta.star.g$lambda     = theta.star.p$lambda[,,perm.indexes[g]]
    theta.star.g$A          = theta.star.p$A[,,perm.indexes[g]]
    for (n in 1:N) {theta.star.g$a[[n]]     = theta.star.p$a[[n]][,perm.indexes[g]]}
    theta.star.g$Gamma      = theta.star.p$Gamma[,,perm.indexes[g]]
    theta.star.g$PR_TR      = theta.star.p$PR_TR[,,perm.indexes[g]]
    theta.star.g$w          = theta.star.p$w[,perm.indexes[g]]
    theta.star.g$xi         = theta.star.p$xi[,,perm.indexes[g]]
    
    
    # step 3.2: ordinate of lambda.star | y, alpha.star 
    qqq.g    = B.SVAR.MSH.reduced.sampler.log.full.lambda(S.aux, priors, restrictions, starting.values=theta.star.g, debug=FALSE)
    qqq.g    = B.SVAR.MSH.short.mcmc(qqq.g,keep.every)
    for (ss in 1:(S.aux/keep.every)){   
      U.SF           = Y %*% t(qqq.g$posteriors$A[,,ss]) - X %*% theta.star.g$Gamma
      lambda.beta.full.tmp    = matrix(NA,N,M-1)
      for(regime in 2:M){
        # Crossproduct of the residuals
        xi.m              = qqq.g$posteriors$S.t[,ss] == regime
        T.m               = sum(xi.m)
        lambda.beta.full  = priors$lambda.beta + T.m
        if (T.m != 0) {
          Sigma.epsilon  = crossprod(U.SF[xi.m,])
          for(equation in 1:N){
            lambda.alpha.full       = priors$lambda.alpha + Sigma.epsilon[equation,equation]
            lambda.beta.full.tmp[equation,regime-1]    = ln.d.inverse.gamma.2(x=theta.star.g$lambda[equation,regime], alpha=lambda.alpha.full, beta=lambda.beta.full)
          }
        }
      }
      log.full.lambda.s[ss,g]       = sum(lambda.beta.full.tmp)
    }
    
#     # step 3.3: ordinate of A.star | y, alpha.star, lambda.star
#     qqq.g    = B.SVAR.MSH.reduced.sampler.log.full.an(S.aux, priors, restrictions, starting.values=theta.star.g, debug=FALSE)
#     qqq.g    = B.SVAR.MSH.short.mcmc(qqq.g,keep.every)
#     st.hist  = MC.hist(qqq.g$posteriors$S.t)
#     theta.star.tmp       = theta.star.g
#     for (ind in 1:max(which(st.hist$count>min.count.no))){
#       theta.star.tmp$xi = diag(M)[,st.hist$unique[,ind]]
#       log.full.A0.s     = c(as.vector(log.full.A0.s),B.SVAR.MSH.A0.log.pdf.vectorized(A0.array = qqq.g$posteriors$A[,,st.hist$which[[ind]]], theta.star = theta.star.tmp, priors, restrictions, S.constant=S.constant))
#     }
#     
#     # step 3.4: ordinate of P.star | y, alpha.star, lambda.star, A.star
#     qqq.g    = B.SVAR.MSH.reduced.sampler.log.full.PR_TR(S.aux, priors, restrictions, starting.values=theta.star.g, debug=FALSE)
#     qqq.g    = B.SVAR.MSH.short.mcmc(qqq.g,keep.every)
#     log.full.PR_TR.s[,g]        = apply(qqq.g$posteriors$S.t,2,B.SVAR.MSH.PR_TR.log.pdf, theta.star.g=theta.star.g, priors=priors, restrictions=restrictions)   
#     t1    = proc.time()
#     if (debug & g%%1==0) {cat(", total time elapsed: ", (t1-t0)[1]/60," [min]")}
  }
  
  lnc.lambda           = quantile(log.full.lambda.s, probs=seq(from=.5, to=.9999999, by=.001))
  log.full.lambda.p    = rep(NA,length(lnc.lambda))
  for (i in 1:length(lnc.lambda)){
    log.full.lambda.p[i]    = lnc.lambda[i] + log(mean(exp(log.full.lambda.s - lnc.lambda[i])))
  }
  log.full.lambda      = log.full.lambda.p[which(is.finite(log.full.lambda.p))[1]]
  
#   log.full.Gamma.s.tmp = apply(log.full.Gamma.s,2,sum)
#   lnc.Gamma            = quantile(log.full.Gamma.s, probs=seq(from=.5, to=.9999999, by=.001))
#   log.full.Gamma.p     = rep(NA,length(lnc.Gamma))
#   for (i in 1:length(lnc.Gamma)){
#     log.full.Gamma.p[i]  = lnc.Gamma[i] + log(mean(exp(log.full.Gamma.s.tmp - lnc.Gamma[i])))
#   }
#   log.full.Gamma       = log.full.Gamma.p[which(is.finite(log.full.Gamma.p))[1]]
#   
#   lnc.A0               = quantile(log.full.A0.s, probs=seq(from=.5, to=.9999999, by=.0001))
#   log.full.A0.p        = rep(NA,length(lnc.A0))
#   for (i in 1:length(lnc.A0)){
#     log.full.A0.p[i]  = lnc.A0[i] + log(mean(exp(log.full.A0.s - lnc.A0[i])))
#   }
#   log.full.A0          = log.full.A0.p[which(is.finite(log.full.A0.p))[1]]
#   
#   log.full.PR_TR       = mean(na.omit(log.full.PR_TR.s)) + log(mean(exp(na.omit(log.full.PR_TR.s) - mean(na.omit(log.full.PR_TR.s)))))
  
  # final computation
#   lnmdd = log.likelihood + log.prior - log.full.Gamma - log.full.lambda - log.full.A0 - log.full.PR_TR
#   cat("\nlog-MDD value: ",lnmdd)
  output         = new.env()
#   output$lnmdd   = lnmdd
  output$components    = list(
#     log.likelihood = log.likelihood,
#     log.prior      = log.prior,
#     log.full.Gamma = log.full.Gamma,
#     log.full.A0    = log.full.A0,
#     log.full.PR_TR = log.full.PR_TR,
    log.full.lambda= log.full.lambda
  )
  output$lnpdf         = list(
#     log.full.Gamma.s  = log.full.Gamma.s.tmp,
#     log.full.A0.s     = log.full.A0.s,
#     log.full.PR_TR.s  = log.full.PR_TR.s,
    log.full.lambda.s = log.full.lambda.s
  )
  
  return(as.list(output))
}










































justforGamma.mdd = function(theta.star, qqq, priors, restrictions, permutation.fraction, S.aux, keep.every=1, debug=TRUE){
   
   cat("\nEstimation of the MDD for the SVAR-MSH model")
   
   S  = dim(qqq$posteriors$lambda)[3]
   N  = dim(qqq$last.draws$Y)[2]
   M  = dim(qqq$last.draws$PR_TR)[1]
   p  = (dim(theta.star$Y)[1] - dim(theta.star$U)[1])
   T  = dim(theta.star$Y)[1]
   TT = T - p
   d     = sum(restrictions$dj) 
   Q     = length(restrictions$dj)
   
   rest  = 0
   for (n in 1:N){
      rest           = rest + sum(diag(restrictions$Q[[n]])==1)
   }
   if (rest!=0){
      restricted     = TRUE
      permutation.fraction    = 1
      cat(" with restricted matrix of contemporaneous effects, A0")
   } else {
      restricted     = FALSE
      cat("\nwith unrestricted matrix of contemporaneous effects, A0")
   }
   
   
   # log of the posterior density ordinate at theta.star
   # step 1: compute permutations
   ##################################################################
   theta.star.p      = B.SVAR.MSH.UNRESTRICTED.state.labels.permutation.augmentation(theta.star)
   if (!restricted){
      theta.star.p      = B.SVAR.MSH.UNRESTRICTED.equation.ordering.permutation.augmentation(theta.star.p)
   }
   permutations.no   = floor(dim(theta.star.p$lambda)[3]*permutation.fraction)
   if (!restricted){
      perm.indexes      = sample.int(n=dim(theta.star.p$lambda)[3], size = permutations.no)
      cat("\nEstimator computed on the basis of ",length(perm.indexes), "random permutations:\n")
   } else {
      perm.indexes      = 1:dim(theta.star.p$lambda)[3]
      cat("\nEstimator computed on the basis of a complete set of ",length(perm.indexes), " permutations:\n")
   }
   
   # step 2: ordinate of alpha.star | y
   ##################################################################
   # X, Y 
   X        = matrix(data=1, nrow=TT, ncol=1)   # A column of ones for the intercept
   if(p > 0){
      for(lag in 1:p){
         Y.m.j    = theta.star$Y[(p+1-lag):(T-lag),]
         X        = cbind(X,Y.m.j)
      }
   }
   Y        = matrix(theta.star$Y[(p+1):T,],ncol=N)   
   
   theta.star.g         = theta.star
   log.full.Gamma.s     = array(NA,c(N,S,length(perm.indexes)))
   
   t0    = proc.time()
   for (g in 1:length(perm.indexes)){
      if (debug & g%%1==0) {cat("\nPermutation No.: ",g)}
      theta.star.g$Sigma      = theta.star.p$Sigma[,,,perm.indexes[g]]
      theta.star.g$lambda     = theta.star.p$lambda[,,perm.indexes[g]]
      theta.star.g$A          = theta.star.p$A[,,perm.indexes[g]]
      for (n in 1:N) {theta.star.g$a[[n]]     = theta.star.p$a[[n]][,perm.indexes[g]]}
      theta.star.g$Gamma      = theta.star.p$Gamma[,,perm.indexes[g]]
      theta.star.g$PR_TR      = theta.star.p$PR_TR[,,perm.indexes[g]]
      theta.star.g$w          = theta.star.p$w[,perm.indexes[g]]
      theta.star.g$xi         = theta.star.p$xi[,,perm.indexes[g]]
      
      # Step 3.1: log-ordinate of alpha.star | y
      # this part uses the posterior draws
      H     = array(NA, dim(priors$Gamma.Var)) 
      P     = vector("list",N)
      label.perm        = allPermutations(M)[perm.indexes[g]%/%prod(1:N) + 1,]
      row.perm          = allPermutations(N)[perm.indexes[g]%%prod(1:N) ,]
      for (s in 1:S){
         # Recreate A0 and lambda after label and row permutations
         A0       = qqq$posteriors$A[,,s]
         la       = qqq$posteriors$lambda[,,s]
         sign.A0  = sign(diag(A0))
         A0.tmp   = diag(sign.A0)%*%A0
         tau1     = 1/diag(A0.tmp)^2
         B0       = diag(sqrt(tau1))%*%A0.tmp
         tau = la
         for (m in 1:M) {
            tau[,m] = la[,m]*tau1
         }
         tau.l    = tau[,label.perm]
         A0.l     = diag(1/sqrt(tau.l[,1]))%*%B0
         A0.l     = diag(sign.A0)%*%A0.l
         lambda.l = tau.l
         for (m in 1:M) {
            lambda.l[,m] = tau.l[,m]/tau.l[,1]
         }
         A0.l     = A0.l[row.perm,]
         lambda.l = lambda.l[row.perm,]
         
         for(equation in 1:N){
            # Crossproducts
            XX   = matrix(0, ncol(X), ncol(X))
            XY   = matrix(0, ncol(X), N)      
            for (regime in 1:M){
               xi.m    = qqq$posteriors$S.t[,s] == regime
               XX      = XX + crossprod(X[xi.m,])/lambda.l[equation,regime]
               XY      = XY + crossprod(X[xi.m,],Y[xi.m,])/lambda.l[equation,regime]
            }
            H[,,equation]     = solve(XX + solve(priors$Gamma.Var[,,equation]))
            P[[equation]]     = H[,,equation] %*% ( XY + solve(priors$Gamma.Var[,,equation]) %*% priors$Gamma.Mean) %*% as.matrix(A0.l[equation,])
            log.full.Gamma.s[equation,s,g] = dmvnorm(x=theta.star.p$Gamma[,equation,g], mean=P[[equation]], sigma=H[,,equation], log=TRUE)
         }  
      }
   }
   
   
   log.full.Gamma.s.tmp = apply(log.full.Gamma.s,2:3,sum)
   log.full.Gamma       = max(log.full.Gamma.s.tmp) + log(mean(exp(log.full.Gamma.s.tmp - max(log.full.Gamma.s.tmp))))
   output         = new.env()
   output$log.full.Gamma = log.full.Gamma
   output$log.full.Gamma.s  = log.full.Gamma.s.tmp
   
   return(as.list(output))
}

















# log.full.lambda   
# lnc = mean(log.full.lambda.s) 
# lnc = max(log.full.lambda.s) 
# lnc = quantile(log.full.lambda.s, probs=.9)
# lnc + log(mean(exp(log.full.lambda.s - lnc)))
# log.full.lambda = lnc + log(mean(exp(log.full.lambda.s - lnc)))
# 
# for ( i in 1:length(sprob)){
#    cat(sprob[i] + log(mean(exp(log.full.lambda.s - sprob[i]))))
# }
# #
# 
# 
# 
# 
# 
# log.full.Gamma       
# lnc = mean(log.full.Gamma.s.tmp)
# lnc = max(log.full.Gamma.s.tmp)
# lnc = quantile(log.full.Gamma.s.tmp, probs=.95)
# lnc
# lnc + log(mean(exp(log.full.Gamma.s.tmp - lnc)))
# #
# 
# 
# 
# 
# 
# 
# S.t    = qqq$posteriors$S.t[,s]
# lambda = qqq$posteriors$lambda[,,s]
# A0     = qqq$posteriors$A[,,s]
# lfgs   = rep(NA,N)
# G      = theta.star.p$Gamma[,,g]
# 
# 
# # S.t    = qqq$posteriors$S.t[,s]
# lambda[1:2,] = lambda[2:1,]
# A0[1:2,]     = A0[2:1,]
# G[,1:2]      = G[,2:1]
# lfgs   = rep(NA,N)
# 
# 
# H     = array(NA, dim(priors$Gamma.Var)) 
# P     = vector("list",N)
# for (s in 1:S){
#    for(equation in 1:N){
#       # Crossproducts
#       XX   = matrix(0, ncol(X), ncol(X))
#       XY   = matrix(0, ncol(X), N)      
#       for (regime in 1:M){
#          xi.m    = S.t == regime
#          XX      = XX + crossprod(X[xi.m,])/lambda[equation,regime]
#          XY      = XY + crossprod(X[xi.m,],Y[xi.m,])/lambda[equation,regime]
#       }
#       H[,,equation]     = solve(XX + solve(priors$Gamma.Var[,,equation]))
#       P[[equation]]     = H[,,equation] %*% ( XY + solve(priors$Gamma.Var[,,equation]) %*% priors$Gamma.Mean) %*% as.matrix(A0[equation,])
#       lfgs[equation] = dmvnorm(x=G[,equation], mean=P[[equation]], sigma=H[,,equation], log=TRUE)
#    }  
# }
# lfgs
# sum(lfgs)















B.SVAR.MSH.short.mcmc = function(qqq,keep.every=10){
   
   if (!is.null(qqq$posteriors$lambda))      S = dim(qqq$posteriors$lambda)[3]
   if (!is.null(qqq$posteriors$A))           S = dim(qqq$posteriors$A)[3]
   if (!is.null(qqq$posteriors$Gamma))       S = dim(qqq$posteriors$Gamma)[3]
   if (!is.null(qqq$posteriors$PR_TR))       S = dim(qqq$posteriors$PR_TR)[3]
   if (!is.null(qqq$posteriors$w))           S = dim(qqq$posteriors$w)[2]
   if (!is.null(qqq$posteriors$S.t))         S = dim(qqq$posteriors$S.t)[2]
   
   N  = dim(qqq$last.draws$Y)[2]
   qqq.new  = qqq
   ind= seq(from=1, to=S, by=keep.every)
   if (!is.null(qqq$posteriors$Sigma))       qqq.new$posteriors$Sigma      = qqq$posteriors$Sigma[,,,ind]
   if (!is.null(qqq$posteriors$lambda))      qqq.new$posteriors$lambda     = qqq$posteriors$lambda[,,ind]
   if (!is.null(qqq$posteriors$A))           qqq.new$posteriors$A          = qqq$posteriors$A[,,ind]
   if (!is.null(qqq$posteriors$a)){
      for (n in 1:N){
         if (dim(qqq$posteriors$a[[n]])[1]==1){
            qqq.new$posteriors$a[[n]]  = matrix(qqq$posteriors$a[[n]][,ind],nrow=1)
         } else {
            qqq.new$posteriors$a[[n]]  = qqq$posteriors$a[[n]][,ind]
         }
      }
   }
   if (!is.null(qqq$posteriors$Gamma))       qqq.new$posteriors$Gamma      = qqq$posteriors$Gamma[,,ind]
   if (!is.null(qqq$posteriors$PR_TR))       qqq.new$posteriors$PR_TR      = qqq$posteriors$PR_TR[,,ind]
   if (!is.null(qqq$posteriors$w))           qqq.new$posteriors$w          = qqq$posteriors$w[,ind]
   if (!is.null(qqq$posteriors$S.t))         qqq.new$posteriors$S.t        = qqq$posteriors$S.t[,ind]
   
   return(qqq.new)
}

B.SVAR.MSH.merge = function(qqq1,qqq2){
   
   N  = dim(qqq1$last.draws$Y)[2]
   
   qqq.new     = qqq1
   qqq.new$posteriors$Sigma      = abind(qqq1$posteriors$Sigma,qqq2$posteriors$Sigma)
   qqq.new$posteriors$lambda     = abind(qqq1$posteriors$lambda,qqq2$posteriors$lambda)
   qqq.new$posteriors$A          = abind(qqq1$posteriors$A,qqq2$posteriors$A)
   
   for (n in 1:N){
      qqq.new$posteriors$a[[n]]      = abind(qqq1$posteriors$a[[n]],qqq2$posteriors$a[[n]])
   }
   qqq.new$posteriors$Gamma      = abind(qqq1$posteriors$Gamma,qqq2$posteriors$Gamma)
   qqq.new$posteriors$PR_TR      = abind(qqq1$posteriors$PR_TR,qqq2$posteriors$PR_TR)
   qqq.new$posteriors$w          = abind(qqq1$posteriors$w,qqq2$posteriors$w)
   qqq.new$posteriors$S.t        = abind(qqq1$posteriors$S.t,qqq2$posteriors$S.t)
   
   return(qqq.new)
}


B.SVAR.MSH.UNRESTRICTED.state.labels.permutation.augmentation = function(aux){
   
   N  = dim(aux$Y)[2]
   M  = dim(aux$PR_TR)[1]
   p  = (dim(theta.star$Gamma)[1] - 1)/N
   
   permutationsM     = allPermutations(M)
   
   qqq.new  = aux
   qqq.new$Sigma   = array(NA,c(N,N,M,prod(1:M)))
   qqq.new$lambda  = array(NA,c(N,M,prod(1:M)))
   qqq.new$A       = array(NA,c(N,N,prod(1:M)))
   for (n in 1:N){
      qqq.new$a[[n]] = array(NA,c(length(aux$a[[n]]),prod(1:M)))
   }
   qqq.new$Gamma   = array(NA,c(1+p*N,N,prod(1:M)))
   qqq.new$PR_TR   = array(NA,c(M,M,prod(1:M)))
   qqq.new$w       = array(NA,c(M*M,prod(1:M)))
   qqq.new$xi      = array(NA,c(dim(aux$xi),prod(1:M)))
   
   A0       = aux$A
   la       = aux$lambda
   
   for (mf in 1:prod(1:M)){
      # Permutation due to label switching
      # create B0 and tau
      
      sign.A0  = sign(diag(A0))
      A0.tmp   = diag(sign.A0)%*%A0
      tau1     = 1/diag(A0.tmp)^2
      B0       = diag(sqrt(tau1))%*%A0   
      tau = la
      for (m in 1:M) {
         tau[,m] = la[,m]*tau1
      }
      
      # Recreate A0 and lambda after label switching
      tau.l    = tau[,permutationsM[mf,]]
      A0.l     = diag(1/sqrt(tau.l[,1]))%*%B0
      A0.l     = diag(sign.A0)%*%A0.l
      lambda.l = tau.l
      for (m in 1:M) {
         lambda.l[,m] = tau.l[,m]/tau.l[,1]
      }
      
      qqq.new$A[,,mf]       = A0.l
      qqq.new$lambda[,,mf]  = lambda.l
      for (n in 1:N){
#          qqq.new$a[[n]][,mf]= A0.l[n,A0.l[n,]!=0]
        qqq.new$a[[n]][,mf]= c(qqq.new$Un[[n]]%*%A0.l[n,])
      }
      qqq.new$Gamma[,,mf]   = aux$Gamma %*% solve(t(A0)) %*% t(A0.l)
      qqq.new$Sigma[,,,mf]  = aux$Sigma[,,permutationsM[mf,]]
      qqq.new$PR_TR[,,mf]   = PR_TR.l = aux$PR_TR[permutationsM[mf,],permutationsM[mf,]]
      qqq.new$w[,mf]        = matrix(t(PR_TR.l),ncol=1)
      qqq.new$xi[,,mf]      = aux$xi[permutationsM[mf,],]
   }
   
   return(qqq.new)
}



B.SVAR.MSH.UNRESTRICTED.equation.ordering.permutation.augmentation = function(qqq.new){
   # This only follows function: B.SVAR.MSH.UNRESTRICTED.state.labels.permutation.augmentation
   N  = dim(qqq.new$Y)[2]
   M  = dim(qqq.new$PR_TR)[1]
   p  = (dim(theta.star$Gamma)[1] - 1)/N
   
   permutationsN     = allPermutations(N)
   qqq.nm  = qqq.new
   qqq.nm$Sigma   = array(NA,c(N,N,M,prod(1:M)*prod(1:N)))
   qqq.nm$lambda  = array(NA,c(N,M,prod(1:M)*prod(1:N)))
   qqq.nm$A       = array(NA,c(N,N,prod(1:M)*prod(1:N)))
   for (n in 1:N){
      qqq.nm$a[[n]] = array(NA,c(dim(qqq.new$a[[n]])[1],prod(1:M)*prod(1:N)))
   }
   qqq.nm$Gamma   = array(NA,c(1+p*N,N,prod(1:M)*prod(1:N)))
   qqq.nm$PR_TR   = array(NA,c(M,M,prod(1:M)*prod(1:N)))
   qqq.nm$w       = array(NA,c(M*M,prod(1:M)*prod(1:N)))
   qqq.nm$xi      = array(NA,c(dim(qqq.new$xi)[1:2],prod(1:M)*prod(1:N)))
   
   
   # Permutation wrt equation ordering
   # Recreate A0 and lambda after equation reordering
   for (s in 1:(prod(1:M))){
      A0       = qqq.new$A[,,s]
      la       = qqq.new$lambda[,,s]
      
      for (mf in 1:prod(1:N)){
         index                               = (s-1)*prod(1:N) + mf
         A0.ro = qqq.nm$A[,,index]       = A0[permutationsN[mf,],]
         qqq.nm$lambda[,,index]  = la[permutationsN[mf,],]
         for (n in 1:N){
            qqq.nm$a[[n]][,index]= A0.ro[n,]
         }
         qqq.nm$Gamma[,,index]   = qqq.new$Gamma[,,s] %*% solve(t(A0)) %*% t(A0.ro)
         qqq.nm$Sigma[,,,index]  = qqq.new$Sigma[,,,s]
         qqq.nm$PR_TR[,,index]   = qqq.new$PR_TR[,,s]
         qqq.nm$w[,index]        = qqq.new$w[,s]
         qqq.nm$xi[,,index]      = qqq.new$xi[,,s]
      }
   }
   return(qqq.nm)
}


B.SVAR.MSH.reduced.sampler.log.full.lambda = function(S, priors, restrictions, starting.values=theta.star.g, debug=FALSE) {
   # reduced Gibbs sampler to be used for the purpose of evaluating the ordinate of the full conditional posterior of lambda
   #----------------------------------------------------------------------------- 
   aux   = starting.values    
   N     = dim(aux$Y)[2]     
   T     = dim(aux$Y)[1]     
   TT    = dim(aux$U)[1]     
   p     = T - TT
   M     = dim(aux$PR_TR)[1]
   d     = sum(restrictions$dj) 
   #-----------------------------------------------------------------------------     
   a     = vector("list")
   for (i in 1:N){
      a[[i]]   = array(NA, c(dim(aux$Un[[i]])[1],S))
   }   
   posteriors = list(        
#       Sigma    = array(NA, c(N,N,M,S)),
#       lambda   = array(NA, c(N,M,S)),
      A        = array(NA, c(N,N,S)),
#       a        = a,
#       Gamma    = array(NA, c(dim(aux$Gamma),S)),
#       PR_TR    = array(NA, c(M,M,S)),
#       w        = array(NA, c(d,S)),
      S.t      = array(NA, c(TT,S))
      #       U        = array(NA, c(dim(aux$U.SF),S)),
      #       order    = array(NA, c(M,S))
   )
   #----------------------------------------------------------------------------- 
   for(iteration in 1:S){
      # Filtering and smoothing step      
      if(debug) print("Filtering and Smoothing step")
      aux = B.SVAR.MSH.filtering.smoothing(aux, iterations=3)
      
      # Hidden Markov Chain step
      if(debug) print("Hidden Markov Chain step") 
      aux = B.SVAR.MSH.hidden.markov.chain(aux, priors, restrictions)
      
      # Inverted-Gamma step
      if(debug) print("Inverted-Gamma step")        
      aux = B.SVAR.MSH.VW.inverted.gamma(aux, priors)
      
      # Structural step with VW algorithm              
      if(debug) print("Structural step")        
      aux = B.SVAR.MSH.VW.structural.VW(aux, priors, restrictions)
      
      # Regression step        
#       if(debug) print("Regression step")        
#       aux = B.SVAR.MSH.regression(aux, priors)
   
      # Save posteriors as vectors        
#       for(regime in 1:M){
#          posteriors$Sigma[,,regime,iteration] = aux$Sigma[,,regime]
#       }
#       posteriors$lambda[,,iteration]   = as.matrix(aux$lambda)
      posteriors$A[,,iteration]        = aux$A
#       for (i in 1:N) {
#          posteriors$a[[i]][,iteration]      = as.matrix(aux$a[[i]])
#       }
#       posteriors$Gamma[,,iteration]    = aux$Gamma
#       posteriors$PR_TR[,,iteration]    = aux$PR_TR
#       posteriors$w[,iteration]         = matrix(aux$w, ncol=1)
      posteriors$S.t[,iteration]       = matrix(max.col(t(aux$xi)), ncol=1)
   }
   #-----------------------------------------------------------------------------    
   output  <- list(                
      last.draws  = aux,                
      posteriors  = posteriors               
   )
   return(output)
}
























B.SVAR.MSH.reduced.sampler.log.full.PR_TR = function(S, priors, restrictions, starting.values=theta.star.g, debug=FALSE) {
   # reduced Gibbs sampler to be used for the purpose of evaluating the ordinate of the full conditional posterior of lambda
   #----------------------------------------------------------------------------- 
   aux   = starting.values    
   N     = dim(aux$Y)[2]     
   T     = dim(aux$Y)[1]     
   TT    = dim(aux$U)[1]     
   p     = T - TT
   M     = dim(aux$PR_TR)[1]
   d     = sum(restrictions$dj) 
   #-----------------------------------------------------------------------------     
   a     = vector("list")
   for (i in 1:N){
      a[[i]]   = array(NA, c(dim(aux$Un[[i]])[1],S))
   }   
   posteriors = list(        
#       Sigma    = array(NA, c(N,N,M,S)),
#       lambda   = array(NA, c(N,M,S)),
#       A        = array(NA, c(N,N,S)),
#       a        = a,
#       Gamma    = array(NA, c(dim(aux$Gamma),S)),
#       PR_TR    = array(NA, c(M,M,S)),
#       w        = array(NA, c(d,S)),
      S.t      = array(NA, c(TT,S))
      #       U        = array(NA, c(dim(aux$U.SF),S)),
      #       order    = array(NA, c(M,S))
   )
   #----------------------------------------------------------------------------- 
   for(iteration in 1:S){
      # Filtering and smoothing step      
      if(debug) print("Filtering and Smoothing step")
      aux = B.SVAR.MSH.filtering.smoothing(aux, iterations=3)
      
      # Hidden Markov Chain step
      if(debug) print("Hidden Markov Chain step") 
      aux = B.SVAR.MSH.hidden.markov.chain(aux, priors, restrictions)
      
#       # Inverted-Gamma step
#       if(debug) print("Inverted-Gamma step")        
#       aux = B.SVAR.MSH.VW.inverted.gamma(aux, priors)
#       
#       # Structural step with VW algorithm              
#       if(debug) print("Structural step")        
#       aux = B.SVAR.MSH.VW.structural.VW(aux, priors, restrictions)
#       
#       # Regression step        
#       if(debug) print("Regression step")        
#       aux = B.SVAR.MSH.regression(aux, priors)
      
      # Save posteriors as vectors        
#       for(regime in 1:M){
#          posteriors$Sigma[,,regime,iteration] = aux$Sigma[,,regime]
#       }
#       posteriors$lambda[,,iteration]   = as.matrix(aux$lambda)
#       posteriors$A[,,iteration]        = aux$A
#       for (i in 1:N) {
#          posteriors$a[[i]][,iteration]      = as.matrix(aux$a[[i]])
#       }
#       posteriors$Gamma[,,iteration]    = aux$Gamma
#       posteriors$PR_TR[,,iteration]    = aux$PR_TR
#       posteriors$w[,iteration]         = matrix(aux$w, ncol=1)
      posteriors$S.t[,iteration]       = matrix(max.col(t(aux$xi)), ncol=1)
   }
   #-----------------------------------------------------------------------------    
   output  <- list(                
      last.draws  = aux,                
      posteriors  = posteriors               
   )
   return(output)
}








B.SVAR.MSH.reduced.sampler.log.full.an = function(S, priors, restrictions, starting.values=theta.star.g, debug=FALSE) {
   # reduced Gibbs sampler to be used for the purpose of evaluating the ordinate of the full conditional posterior of lambda
   #----------------------------------------------------------------------------- 
   aux   = starting.values    
   N     = dim(aux$Y)[2]     
   T     = dim(aux$Y)[1]     
   TT    = dim(aux$U)[1]     
   p     = T - TT
   M     = dim(aux$PR_TR)[1]
   d     = sum(restrictions$dj) 
   #-----------------------------------------------------------------------------     
   a     = vector("list")
   for (i in 1:N){
      a[[i]]   = array(NA, c(dim(aux$Un[[i]])[1],S))
   }   
   posteriors = list(        
      #       Sigma    = array(NA, c(N,N,M,S)),
      #       lambda   = array(NA, c(N,M,S)),
      A        = array(NA, c(N,N,S)),
      a        = a,
      #       Gamma    = array(NA, c(dim(aux$Gamma),S)),
      #       PR_TR    = array(NA, c(M,M,S)),
      #       w        = array(NA, c(d,S)),
      S.t      = array(NA, c(TT,S))
      #       U        = array(NA, c(dim(aux$U.SF),S)),
      #       order    = array(NA, c(M,S))
   )
   #----------------------------------------------------------------------------- 
   for(iteration in 1:S){
      # Filtering and smoothing step      
      if(debug) print("Filtering and Smoothing step")
      aux = B.SVAR.MSH.filtering.smoothing(aux, iterations=3)
      
      # Hidden Markov Chain step
      if(debug) print("Hidden Markov Chain step") 
      aux = B.SVAR.MSH.hidden.markov.chain(aux, priors, restrictions)
      
      #             # Inverted-Gamma step
      #             if(debug) print("Inverted-Gamma step")        
      #             aux = B.SVAR.MSH.VW.inverted.gamma(aux, priors)
      
      # Structural step with VW algorithm              
      if(debug) print("Structural step")        
      aux = B.SVAR.MSH.VW.structural.VW(aux, priors, restrictions)
      
      #             # Regression step        
      #             if(debug) print("Regression step")        
      #             aux = B.SVAR.MSH.regression(aux, priors)
      
      # Save posteriors as vectors        
      #       for(regime in 1:M){
      #          posteriors$Sigma[,,regime,iteration] = aux$Sigma[,,regime]
      #       }
      #       posteriors$lambda[,,iteration]   = as.matrix(aux$lambda)
      posteriors$A[,,iteration]        = aux$A
      for (i in 1:N) {
         posteriors$a[[i]][,iteration]      = as.matrix(aux$a[[i]])
      }
      #       posteriors$Gamma[,,iteration]    = aux$Gamma
      #       posteriors$PR_TR[,,iteration]    = aux$PR_TR
      #       posteriors$w[,iteration]         = matrix(aux$w, ncol=1)
      posteriors$S.t[,iteration]       = matrix(max.col(t(aux$xi)), ncol=1)
   }
   #-----------------------------------------------------------------------------    
   output  <- list(                
      last.draws  = aux,                
      posteriors  = posteriors               
   )
   return(output)
}









###############################################################
# NOT NEEDED
###############################################################
# B.SVAR.MSH.VW.structural.VW.reduced = function(n.reduced, aux, priors, restrictions){
#    # Setup constants     
#    #-----------------------------------------------------------------------------    
#    N  = dim(aux$Y)[2]     
#    T  = dim(aux$Y)[1]     
#    TT = dim(aux$U)[1]     
#    p  = T - TT
#    # X, Y 
#    #-----------------------------------------------------------------------------
#    X     = matrix(data=1, nrow=TT, ncol=1)   # A column of ones for the intercept
#    if(p > 0){
#       for(lag in 1:p){
#          Y.m.j    = aux$Y[(p+1-lag):(T-lag),]
#          X        = cbind(X,Y.m.j)
#       }
#    }
#    Y     = matrix(aux$Y[(p+1):T,],ncol=N)      
#    X     = t(X)
#    Y     = t(Y)
#    AA.sm1= aux$A # The draw of A used later for the normalization
#    # Construct Omega_a_n and mu_a_n
#    #----------------------------------------------------------------------------- 
#    Omega.bar      = vector("list", N)
#    mu.bar         = vector("list", N)
#    Omega.bar.tmp  = vector("list", N)
#    mu.bar.tmp     = vector("list", N)
#    
#    for (n in 1:N){
#       Omega.bar.tmp[[n]]   = matrix(0, dim(aux$Un[[n]])[1], dim(aux$Un[[n]])[1])
#       mu.bar.tmp[[n]]      = matrix(0, 1, dim(aux$Un[[n]])[1])
#       for (regime in 1:M){
#          xi.m                 = aux$xi[regime,]==1
#          Omega.bar.tmp[[n]]   = Omega.bar.tmp[[n]] + tcrossprod(aux$Un[[n]]%*%Y[,xi.m]/sqrt(aux$lambda[n,regime]))
#          mu.bar.tmp[[n]]      = mu.bar.tmp[[n]]    + (aux$Gamma[,n] %*% X[,xi.m] %*% t(aux$Un[[n]] %*% Y[,xi.m])/aux$lambda[n,regime])
#       }
#       Omega.bar[[n]] = solve( Omega.bar.tmp[[n]] + solve(aux$Un[[n]] %*% priors$A.Var[[n]] %*% t(aux$Un[[n]])) + aux$Un[[n]] %*% t(priors$Gamma.Mean) %*% solve(priors$Gamma.Var[,,n]) %*% priors$Gamma.Mean %*% t(aux$Un[[n]]))
#       mu.bar[[n]]    = (mu.bar.tmp[[n]] + aux$Gamma[,n]%*%solve(priors$Gamma.Var[,,n])%*%priors$Gamma.Mean%*%t(aux$Un[[n]]) ) %*% Omega.bar[[n]]
#    }
#    # Algorithm from Waggoner Zha (2003): RWZrestrictions_FiscalPolicy \ fn_gibbsrvar.m
#    #----------------------------------------------------------------------------- 
#    w     = array(0, c(N,1))
#    
#    for(i.star in n.reduced:N){     # given last A0gbs and generate new A0bgs  
#       A.star   = t(aux$A[-i.star,])
#       w        = orthogonal.complement.matrix.TW(A.star)
#       #*** Constructing orthonormal basis w_1, ..., w_qi at each Gibbs step
#       Tn    = t(chol(TT*Omega.bar[[i.star]]))    # defined below eq (14)
#       w0    = t(Tn) %*% aux$Un[[i.star]] %*% w
#       w1    = w0 / as.numeric(sqrt(crossprod(w0)))
#       if (dim(aux$Un[[i.star]])[1]>1) {
#          W     = orthogonal.complement.matrix.TW(w1)
#       }
#       # compute the means of the Absolute Normal and normal distributions distribution:
#       xi    = rep(NA,dim(aux$Un[[i.star]])[1])
#       xi[1] = mu.bar[[i.star]] %*% solve(t(Tn)) %*% w1
#       if (length(xi)>1){
#          for (j in 2:length(xi)){
#             xi[j]    = mu.bar[[i.star]] %*% solve(t(Tn)) %*% W[,j-1]
#          }
#       }
#       #*** Draw beta's according to Proposition C.1 of Villani(2009, JAE).
#       # All but first element
#       gkbeta         = array(NA, c(dim(aux$Un[[i.star]])[1],1))  # qi-by-1: greak beta's
#       if (length(xi)>1){
#          for (j in 2:length(xi)){
#             gkbeta[j,]    = rnorm(n=1, mean=xi[j], sd=sqrt(1/TT))
#          }
#       }
#       # Sample first element beta1
#       # A draw from the absolute normal distribution according to 
#       mu.an       = rep(NA,2)
#       s.an        = rep(NA,2)
#       mu.an[1]    = .5*xi[1] - .5*sqrt(xi[1]^2+4)
#       mu.an[2]    = .5*xi[1] + .5*sqrt(xi[1]^2+4)
#       for (j in 1:2) s.an[j]  = ((mu.an[j]^2)*(1/TT))/(1+mu.an[j]^2)
#       www         = 1/(1+exp(2*xi[1]*TT))
#       param       = norMix(mu=mu.an, sigma = sqrt(s.an), w = c(www,1-www))
#       gkbeta[1,]  = rnorMix(n=1, obj=param)
#       
#       if (length(xi)>1) {
#          WWW = cbind(w1,W)
#       } else {
#          WWW = w1
#       }
#       aux$a[[i.star]]   = as.vector(Tn %*% WWW %*% gkbeta)    # equation (15)
#       aux$A[i.star,]    = aux$a[[i.star]] %*% aux$Un[[i.star]]
#    }
#    inv.A          = solve(aux$A)   
#    for (regime in 1:M){
#       aux$Sigma[,,regime]     = inv.A %*% diag(aux$lambda[,regime]) %*% t(inv.A)
#    }
#    return(aux)
# }






























# B.SVAR.MSH.permutation = function(aux, permutation="random", which.lambda=1){
#    # argument: permutation = c("random", "ordered", "none")
#    # setup constants
#    N  = dim(aux$Y)[2]     
#    T  = dim(aux$Y)[1]     
#    TT = dim(aux$U)[1]     
#    p  = T - TT 
#    d  = sum(restrictions$dj) 
#    Q  = length(restrictions$dj)
#    
#    if (which.lambda>M) which.lambda = 1
#    
#    if (permutation=="none") {} else {
#       if (permutation=="random") {
#          new.perm       = shuffle(M)
#       } else if (permutation=="ordered") {
#          new.perm       = order(aux$lambda[which.lambda,])
#       }
#       
#       # premute
#       aux$Sigma   = aux$Sigma[,,new.perm]
#       aux$lambda  = aux$lambda[,new.perm]
#       aux$PR_TR   = aux$PR_TR[new.perm,new.perm]
#       aux$xi      = aux$xi[new.perm,]            
#    } 
#    
#    return(aux)
# }