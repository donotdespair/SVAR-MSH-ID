#---------------------------------------------------------------------------------------------------
#  B.VECM: Useful functions
#---------------------------------------------------------------------------------------------------


library(permute)
library(abind)

dig2     = function(x,a,b,logarithm=FALSE){
   lcig2    = lgamma(a/2) + (-a/2)*log(b/2)
   lkig2    = (-(a+2)/2)*log(x) - 0.5*(b/x)
   out      = as.numeric(lkig2 - lcig2)
   if (!logarithm) out = exp(out)
   return(out)
}

dt2      = function(x,mu,precision,a,b,logarithm=FALSE){
   N        = length(x)
   lct2     = lgamma(a/2) - lgamma((a+N)/2) + (N/2)*log(pi) - 0.5*determinant(precision, logarithm=TRUE)$modulus - (a/2)*log(b)
   lkt2     = (-(a+N)/2)*log(b + t(x - mu) %*% precision %*% (x - mu))
   out      = as.numeric(lkt2 - lct2)
   if (!logarithm) out = exp(out)
   return(out)
}





allPermutations = function(n){
   return(rbind(1:n,allPerms(n)))
}

orthogonal.complement.matrix.TW = function(x){
   # x is a mxn matrix and m>n
   # the function returns a mx(m-n) matrix, out, that is an orthogonal complement of x, i.e.:
   # t(x)%*%out = 0 and det(cbind(x,out))!=0
   N     = dim(x)
   tmp   = qr.Q(qr(x, tol = 1e-10),complete=TRUE)
   out   = as.matrix(tmp[,(N[2]+1):N[1]])
   return(out)
}

B.VECM.Square.Root = function(A, inverse=F){

   # The square root of matrices is computed as in Abadir, Magnus (2005) p. 220-221
   A.eigen = eigen(A)
   if (dim(A)[1] == 1) {
      Lambda = as.vector(sqrt(A.eigen$values))
   } else {
      Lambda = diag(as.vector(sqrt(A.eigen$values)))
   }
   
   if (inverse){
      Lambda = solve(Lambda)
   }
   output  = A.eigen$vectors %*% Lambda %*% t(A.eigen$vectors)
   
   # Output     
   return(output)
}


Null.TW = function(x){
   ww    = as.matrix(abs(Null(x)))
   index = order(apply(ww,2,function(y){which(y==1)}))
   return(as.matrix(ww[,index]))
}

#---------------------------------------------------------------------------------------------------
# Gibbs: Probabilities of regime from posteriors
#---------------------------------------------------------------------------------------------------
# posteriors: from the Gibbs sampler
# burning: disregard the first estimations of the Gibbs

# G.likelihood.function = function(aux){
#    # Computation of the likelihood function for MS-VARs is presented in 
#    # Fruhwirth-Schnatter (2006, book): Section 11.4.1: second eq. on p.332 and eq. (11.3)
#    #-----------------------------------------------------------------------------------------------
#       eta.t    <- matrix(NA, M, TT)
#       # Residuals calculated from Beta.star: depending only on Beta.star, they can be computed now
#       for(regime in 1:M){
#          # Residuals: from Beta star, used for additional Gibbs runs
#          aux$Beta[,,regime]  <- matrix(posteriors$Beta[,regime,iteration], ncol=K)
#          aux$Sigma[,,regime] <- matrix(posteriors$Sigma[,regime,iteration], ncol=K)
#          res.tmp             <- Y - X %*% aux$Beta[,,regime]
#          aux$U[,,regime]     <- matrix(res.tmp, ncol=K)   
#          
#          # Densities conditional on regimes: p(y_t|s_t=m,\mathbf{y}_{t-1},\theta)
#          eta.tmp         <- dmvnorm(matrix(aux$U[,,regime],ncol=K), sigma=matrix(aux$Sigma[,,regime],ncol=K)) 
#          eta.t[regime,]  <- matrix(eta.tmp, nrow=1)
#          
#          # Filtered probabilities: Pr[s_t=m|\mathbf{y}_{t-1},\theta]
#          ## Rebuild drawn regime history
#          aux$xi[regime,] <- 0
#          aux$xi[regime,which(posteriors$S.t[,iteration]==regime)] <- 1
#       }
#       
#       PR_TR.it     = matrix(posteriors$PR_TR[,iteration],ncol=M)
#       PR.filtered  = BHLK.filter(u=aux$U, PR_TR=PR_TR.it, Sigma=aux$Sigma, xi.bar=Ergodic.PR_TR(PR_TR.it))
#       
#       # Log likelihood 
#       log.likelihood[iteration]   <- sum( log(colSums(eta.t * PR.filtered$fp)))
# }

G.probabilities.regime <- function(posteriors, keep.iterations){
   
   # Setup constants
   S.t     <- posteriors$S.t[,keep.iterations]
   T       <- dim(S.t)[1]
   N       <- dim(S.t)[2]
   M       <- sqrt(dim(posteriors$PR_TR)[1])
   
   # Matrix of probabilities
   P.S.t   <- array(NA,c(T,M))
   
   # Function for counting occurrences of regimes    
   count.regimes <- function(S.t, regime){
      output <- sum(S.t[which(S.t == regime)]) / regime
      return(output)
   }
   
   # Using apply()
   for(regime in 1:M){
      P.S.t[,regime]  <- apply(S.t, 1, count.regimes, regime=regime) / N
   }                                             
   
   # Output
   return(P.S.t)
}



#---------------------------------------------------------------------------------------------------
# Plot generated vs. estimated probabilities
#---------------------------------------------------------------------------------------------------
G.plot.probabilities  <- function(Gibbs.output, Xi, keep.iterations){
   
   # Compute estimated probabilities
   Est.P   <- G.probabilities.regime(Gibbs.output$posteriors,keep.iterations)
   
   # Setup constants     
   TT      <- dim(Est.P)[1] 
   M       <- dim(Est.P)[2]
   p       <- length(Xi)[1] - TT
   
   # Remove first p observations from the generated vector of regimes 
   Xi      <- Xi[-(1:p)] 
   
   # Transform the vector of generated regimes into a matrix of probabilities (T x M)
   Gen.P   <- array(NA, c(TT,M))
   
   # Function for counting occurrences of regimes    
   vector.regime <- function(Xi, regime){
      if(Xi==regime){
         output  <- 1
      }else{
         output  <- 0
      }
      return(output)
   }
   
   # Using apply()
   for(regime in 1:M){
      Gen.P[,regime]  <- apply(as.matrix(Xi), 1, vector.regime, regime=regime)
   }
   
   
   # plot smoothed probabilities against known Xi
   op <- par(mfrow = c(M,1), mar=c(2,1,3,0.5), oma=c(1.5,2,1,1)) # M graphs + margins settings
   
   for(regime in 1:M){
      series.plot <- ts(cbind(Est.P[,regime], Gen.P[,regime]))
      title   <- paste("Drawn probabilities, regime ", regime, sep="")
      plot(series.plot, main=title ,ylim=c(0, 1), plot.type="single",col=c(1,2),lty=c(1,2))  
   }
   par(op) # At end of plotting, reset to previous settings
}                               




B.SVAR.MSH.identified.short.mcmc = function(qqq,keep.every=10){
   
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
   if (!is.null(qqq$posteriors$lambda.omega))      qqq.new$posteriors$lambda.omega     = qqq$posteriors$lambda.omega[,,ind]
   if (!is.null(qqq$posteriors$A))           qqq.new$posteriors$A          = qqq$posteriors$A[,,ind]
   if (!is.null(qqq$posteriors$a))           qqq.new$posteriors$a          = qqq$posteriors$a[,ind]
   if (!is.null(qqq$posteriors$alpha))       qqq.new$posteriors$alpha      = qqq$posteriors$alpha[,,ind]
   if (!is.null(qqq$posteriors$PR_TR))       qqq.new$posteriors$PR_TR      = qqq$posteriors$PR_TR[,,ind]
   if (!is.null(qqq$posteriors$hyper))           qqq.new$posteriors$hyper          = qqq$posteriors$hyper[,ind]
   if (!is.null(qqq$posteriors$w))           qqq.new$posteriors$w          = qqq$posteriors$w[,ind]
   if (!is.null(qqq$posteriors$S.t))         qqq.new$posteriors$S.t        = qqq$posteriors$S.t[,ind]
   if (!is.null(qqq$posteriors$U.SF))       qqq.new$posteriors$U.SF      = qqq$posteriors$U.SF[,,ind]
   if (!is.null(qqq$posteriors$MH.logkernel))       qqq.new$posteriors$MH.logkernel      = qqq$posteriors$MH.logkernel[ind]
   
   return(qqq.new)
}

B.SVAR.MSH.identified.merge = function(qqq1,qqq2){
   
   N  = dim(qqq1$last.draws$Y)[2]
   
   qqq.new     = qqq1
   qqq.new$posteriors$Sigma      = abind(qqq1$posteriors$Sigma,qqq2$posteriors$Sigma)
   qqq.new$posteriors$lambda.omega     = abind(qqq1$posteriors$lambda.omega,qqq2$posteriors$lambda.omega)
   qqq.new$posteriors$A          = abind(qqq1$posteriors$A,qqq2$posteriors$A)
   qqq.new$posteriors$a          = abind(qqq1$posteriors$a,qqq2$posteriors$a)
   qqq.new$posteriors$alpha      = abind(qqq1$posteriors$alpha,qqq2$posteriors$alpha)
   qqq.new$posteriors$PR_TR      = abind(qqq1$posteriors$PR_TR,qqq2$posteriors$PR_TR)
   qqq.new$posteriors$hyper      = abind(qqq1$posteriors$hyper,qqq2$posteriors$hyper)
   qqq.new$posteriors$w          = abind(qqq1$posteriors$w,qqq2$posteriors$w)
   qqq.new$posteriors$S.t        = abind(qqq1$posteriors$S.t,qqq2$posteriors$S.t)
   qqq.new$posteriors$U.SF        = abind(qqq1$posteriors$U.SF,qqq2$posteriors$U.SF)
   qqq.new$posteriors$MH.logkernel        = c(qqq1$posteriors$MH.logkernel,qqq2$posteriors$MH.logkernel)
   return(qqq.new)
}
