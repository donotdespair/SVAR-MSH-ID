
B.SVAR.identified.MDD = function(Sp=1000,Gibbs.output, priors, restrictions){ #, debug=FALSE){

    #----------------------------------------------------------------------------- 
    # Setup constants 
    #-----------------------------------------------------------------------------
    aux         <- Gibbs.output$last.draws
    posteriors  <- Gibbs.output$posteriors
    p           <- (dim(aux$Y)[1] - dim(aux$U)[1])
    T           <- dim(aux$Y)[1]
    TT          <- T - p
    N           <- dim(aux$Y)[2] 
    S           <- dim(posteriors$alpha)[3]
    r           = dim(restrictions$Q)[2]

    # Y and X matrices
    XX   <- matrix(data=1, nrow=TT, ncol=1)	# A column of ones for the intercept
    if(p > 0){
        for(lag in 1:p){
            Y.m.j   <- aux$Y[(p+1-lag):(T-lag),]
            XX       <- cbind(XX,Y.m.j)
        }
    }
    Y <- matrix(aux$Y[(p+1):T,],ncol=N)  

    # form a matrix of posterior draws:
    posteriormatrix    = matrix(NA, S, 0)
    posteriormatrix       = cbind(posteriormatrix,t(posteriors$lambda.omega))
    
    posteriormatrix    = cbind(posteriormatrix,t(posteriors$a))
    for (n in 1:N){
       posteriormatrix       = cbind(posteriormatrix,t(posteriors$alpha[,n,]))
    }
    posteriormatrix    = as.matrix(cbind(posteriormatrix,t(posteriors$hyper)))
    
    # Likelihood
    log.likelihood   = apply(X=posteriormatrix, MARGIN=1, FUN=loglikelihoodfordraws, restrictions=restrictions,N=N,r=r,Y=Y,XX=XX) # log-likelihood at draws
    post.mean        = apply(posteriormatrix,2,mean)
    post.cov         = cov(posteriormatrix)
    post.cov         = 0.5*(post.cov + t(post.cov))
    lowerb           = c(rep(0,N), rep(-Inf,r + N*(1+N*p)), rep(0, 3) )
    upperb           = c(rep(Inf,N + r + N*(1+N*p)), rep(Inf,3) )
       
    is.draws         = rmvnorm(n=Sp, mean = post.mean,  sigma = post.cov)
    id.in.range      = apply(X=is.draws, MARGIN=1, FUN=in.range, lb=lowerb, ub=upperb)
    is.draws         = is.draws[id.in.range,]

    ll               = apply(X=is.draws, MARGIN=1, FUN=loglikelihoodfordraws,restrictions=restrictions,N=N,r=r,Y=Y,XX=XX) # log-likelihood at draws
    Indicator.setA   = (ll >= min(log.likelihood))
    ls               = apply(is.draws,1,dmvnorm, mean=post.mean, sigma=post.cov, log=TRUE) # log-importance density at draws
    lp               = apply(is.draws,1,log.prior.for.is.draws) # log-prior density at draws
    lis              = ll[Indicator.setA] + lp[Indicator.setA] - ls[Indicator.setA]
    
    mlis             = max(lis)
    log.mdd          = mlis + log(sum(exp(lis - mlis))) - log(length(ll))
    
    cat("Fraction of IS draws in set A: ",sum(Indicator.setA)/Sp)
    
    output           = new.env()
    output$log.mdd   = log.mdd
    output$log.kernel.atIS       = lis
    output$log.likelihood.atIS   = ll
    output$log.prior.atIS        = lp
    output$log.is.atIS           = ls
    output$indicator.setA        = Indicator.setA
    return(as.list(output))
} 
    



#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

loglikelihoodfordraws = function(draws,restrictions,N,r,Y,XX){
   lambda.omega     = matrix(draws[1:N], ncol=1)
   a                = matrix(draws[(N + 1):(N + r)], ncol=1)
   A                = matrix(restrictions$Q %*% a + restrictions$q, ncol=N)
   A.inv            = solve(A)
   Sigma            = array(NA, c(N,N))
   lambda           = as.vector(lambda.omega[,1])
   Sigma            = A.inv %*% diag(lambda) %*% t(A.inv)
   
   # cat("\na")
   alpha            = matrix(draws[(N + r + 1):(N + r + N*(1+N*p))], ncol=N)
   hyper            = matrix(draws[(N + r + N*(1+N*p) + 1):(N + r + N*(1+N*p) + 3)], ncol=1)
   # cat("\nb")
   U.SF             = Y %*% t(A) - XX %*% alpha
   U                = U.SF %*% t(A.inv)
   # cat("\nc")
   # Densities for each regimes
   eta.t 	<- matrix(NA, 1, TT)
   eta.tmp         <- dmvnorm(U, sigma=Sigma) 

   # cat("\nd")
   aux.tmp          = new.env()
   aux.tmp$Y        = matrix(NA, T, N)
   aux.tmp$Sigma    = Sigma
   aux.tmp$U        = U
   aux.tmp          = as.list(aux.tmp)

   output           = as.numeric(sum( log(eta.tmp)))
   
   return(output)
}

log.prior.for.is.draws = function(draws){
   
   lambda.omega     = matrix(draws[1:N], ncol=1)
   a                = matrix(draws[(N + 1):(N + r)], ncol=1)
   A                = matrix(restrictions$Q %*% a + restrictions$q, ncol=N)
   alpha            = matrix(draws[(N + r + 1):(N + r + N*(1+N*p))], ncol=N)
   hyper            = matrix(draws[(N + r + N*(1+N*p) + 1):(N + r + N*(1+N*p) + 3)], ncol=1)
   
   lp.mu          = dt2(x=alpha[1,], mu=rep(0,N), precision=diag(N), a=priors$hyper.a, b=priors$hyper.b, logarithm=TRUE)
   lp.beta        = 0
   for (n in 1:N){
      lp.beta     = lp.beta + dt2(x=alpha[-1,n], mu=as.vector(A[n,]%*%t(priors$beta.P)), precision=priors$beta.H, a=priors$hyper.a, b=priors$hyper.b, logarithm=TRUE)
   }
   lp.a           = dt2(x=as.vector(a), mu=rep(0,r), precision=diag(r), a=priors$hyper.a, b=priors$hyper.b, logarithm=TRUE)
   lp.lambda      = sum(dig2(x=lambda.omega[,1], a=priors$lambda.a, b=priors$lambda.b, logarithm=TRUE))
   lp.hyper       = sum(dig2(x=as.vector(hyper), a=priors$hyper.a, b=priors$hyper.b, logarithm=TRUE))

   lp             = lp.mu + lp.beta + lp.a + lp.lambda + lp.hyper
   return(as.numeric(lp))
}


in.range = function(draws,lb,ub){
   q = as.logical(prod(draws > lb & draws < ub))
   return(q)
}













#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------



MSVAR.IRF.posteriors.matrix <- function(Gibbs.output, restrictions)
{
   #----------------------------------------------------------------------------- 
   # Setup constants 
   #-----------------------------------------------------------------------------
   aux         <- Gibbs.output$last.draws
   posteriors  <- Gibbs.output$posteriors
   M           <- dim(aux$PR_TR)[1]
   p           <- (dim(aux$Y)[1] - dim(aux$U)[1])
   T           <- dim(aux$Y)[1]
   TT          <- T - p
   N           <- dim(aux$Y)[2] 
   S           <- dim(posteriors$Gamma)[3]
   Q           <- length(restrictions$dj)    
   
   # Vectorize all the posteriors for all Gibbs iterations: number of parameters * N
   # PR_TR
   number.parameters   <- sum(restrictions$dj -1)
   cat("\nPR_TR:", number.parameters, "\n")
   # lambda
   number.parameters   <- number.parameters + (M-1) * N
   cat("lambda:", number.parameters, "\n")
   # A
   A.restrict = 0
   for (i in 1:N){
      A.restrict = A.restrict + length(which(diag(restrictions$Q[[i]])==1))
   }
   number.parameters   <- number.parameters + N^2 - A.restrict
   cat("A:", number.parameters, "\n")
   # Gamma
   number.parameters   <- number.parameters + N+p*(N^2)
   cat("Gamma:", number.parameters, "\n")
   
   # Posterior vectors
   posteriors.vector   <- array(NA, c(0,S))
   
   # PR_TR
   for(j in 1:Q){
      index.begin <- sum(restrictions$dj[1:j-1]) + 1
      dj          <- restrictions$dj[j]
      index.end   <- index.begin + dj - 1
      posteriors.vector <- rbind(posteriors.vector,posteriors$w[index.begin:(index.end-1),])
   }  
   
   # lambda
   for (regime in 2:M){
      posteriors.vector = rbind(posteriors.vector, posteriors$lambda[,regime,])
   }
   
   # A
   for (equation in 1:N){
      posteriors.vector = rbind(posteriors.vector, posteriors$a[[equation]])
   }
   
   # Gamma
   for (i in 1:N){
      posteriors.vector = rbind(posteriors.vector, posteriors$Gamma[,i,])
   }
   
   
   #-----------------------------------------------------------------------------------------------
   # output
   #-----------------------------------------------------------------------------------------------
   output  <- list(    posteriors.vector = posteriors.vector,
                       number.parameters = number.parameters   )
   
   return(output)      
   
}




MSVAR.IRF.remove.posteriors.frequency   <- function(Gibbs.output, keep.every)
{
    # parameters
    #-----------------------------------------------------------------------------
    aux <- Gibbs.output$last.draws
    M   <- dim(aux$PR_TR)[1]
    p   <- (dim(aux$Y)[1] - dim(aux$U)[1])
    T   <- dim(aux$Y)[1]
    TT  <- T - p
    K   <- dim(aux$Y)[2]
    d   <- sum(restrictions$dj)   
    N   <- dim(Gibbs.output$posteriors$Gamma)[3] 

    # List of vectorized posteriors for each iteration
    #----------------------------------------------------------------------------- 
    new.N   <- trunc(N/keep.every)
    a          = vector("list",K)
    for (i in 1:K){
       a[[i]]  = array(NA, c(dim(Gibbs.output$posteriors$a[[i]])[1],new.N))
    }
    posteriors <- list(
       Sigma   = array(NA,c(K,K,M,new.N)),
       lambda  = array(NA,c(K,M,new.N)),
       A       = array(NA,c(K,K,new.N)),
       a       = a,
       Gamma   = array(NA,c(1+p*K,K,new.N)),
       PR_TR   = array(NA,c(M,M,new.N)),
       w      = array(NA, c(d,new.N)),
       S.t     = array(NA, c(TT,new.N))
    )
       
    # Build new posterior structure
    #-----------------------------------------------------------------------------
    for(iteration in 1:new.N){
        ### print(iteration*keep.every)
        posteriors$Sigma[,,,iteration]    <- Gibbs.output$posteriors$Sigma[,,,iteration*keep.every]
        posteriors$lambda[,,iteration]   <- Gibbs.output$posteriors$lambda[,,iteration*keep.every]
        posteriors$A[,,iteration]    <- Gibbs.output$posteriors$A[,,iteration*keep.every]
        for (i in 1:K){
           posteriors$a[[i]][,iteration]    <- Gibbs.output$posteriors$a[[i]][,iteration*keep.every]
        }
        posteriors$Gamma[,,iteration]    <- Gibbs.output$posteriors$Gamma[,,iteration*keep.every]
        posteriors$PR_TR[,,iteration]    <- Gibbs.output$posteriors$PR_TR[,,iteration*keep.every]
        posteriors$w[,iteration]        <- Gibbs.output$posteriors$w[,iteration*keep.every]
        posteriors$S.t[,iteration]      <- Gibbs.output$posteriors$S.t[,iteration*keep.every]
    
    }

    # Output
    #-----------------------------------------------------------------------------
    output  <- list(
                last.draws  = aux,
                posteriors  = posteriors
               )

    return(output)  
}



#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------




MSVAR.IRF.remove.posteriors.burnin   <- function(Gibbs.output, burnin.length)
{
    # parameters
    #-----------------------------------------------------------------------------
    aux <- Gibbs.output$last.draws
    M   <- dim(aux$PR_TR)[1]
    p   <- (dim(aux$Y)[1] - dim(aux$U)[1])
    T   <- dim(aux$Y)[1]
    TT  <- T - p
    K   <- dim(aux$Y)[2]
    d   <- sum(restrictions$dj)   
    N   <- dim(Gibbs.output$posteriors$Beta)[3] 

    # List of vectorized posteriors for each iteration
    #----------------------------------------------------------------------------- 
    new.N   <- N - burnin.length
    posteriors <- list(
        Beta    = array(NA, c(K*(1+K*p),M,new.N)),
        Sigma   = array(NA, c(K*K,M,new.N)),
        PR_TR   = array(NA, c(M*M,new.N)),
        w   	= array(NA, c(d,new.N)),
        S.t     = array(NA, c(TT,new.N)),
        F.S.t   = array(NA, c(h,new.N)),
        F.Y     = array(NA, c(h,K,new.N))
        )

    # Build new posterior structure
    #-----------------------------------------------------------------------------
    for(iteration in 1:new.N){
        ### print(iteration*keep.every)
        posteriors$Beta[,,iteration]    <- Gibbs.output$posteriors$Beta[,,iteration+burnin.length]
        posteriors$Sigma[,,iteration]   <- Gibbs.output$posteriors$Sigma[,,iteration+burnin.length]
        posteriors$PR_TR[,iteration]    <- Gibbs.output$posteriors$PR_TR[,iteration+burnin.length]
        posteriors$w[,iteration]        <- Gibbs.output$posteriors$w[,iteration+burnin.length]
        posteriors$S.t[,iteration]      <- Gibbs.output$posteriors$S.t[,iteration+burnin.length]
        posteriors$F.S.t[,iteration]    <- Gibbs.output$posteriors$F.S.t[,iteration+burnin.length]
        posteriors$F.Y[,,iteration]     <- Gibbs.output$posteriors$F.Y[,,iteration+burnin.length]
    
    }

    # Output
    #-----------------------------------------------------------------------------
    output  <- list(
                last.draws  = aux,
                posteriors  = posteriors
               )

    return(output)  
}





 
#---------------------------------------------------------------------------------------------------
# diwish calls encounter overflows here, due to a high degrees of freedom
# logarithm of the function, adapted from MCMCpack, trick inspired from the MSBVAR package
#---------------------------------------------------------------------------------------------------
ldiwish <- function (W, v, S) 
{
    k <- nrow(S)
    
    # denominator 
    lgammapart <- 0
    for(i in 1:k) {
        lgammapart  <- lgammapart + lgamma((v + 1 - i)/2)
    }
    denom <- lgammapart + log(2)*(v * k/2) + log(pi) * (k * (k - 1)/4)
  
    # numerator
    detS <- determinant(S)$modulus[1] 
    detW <- determinant(W)$modulus[1]  
    hold <- S %*% solve(W)
    tracehold <- sum(hold[row(hold) == col(hold)])
    num <- detS*(v/2) + detW*(-(v + k + 1)/2) + (-1/2 * tracehold)

    return(num - denom)
}          
