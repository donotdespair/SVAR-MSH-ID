B.SVAR.MSH.identified.Gibbs = function(S, priors, restrictions, starting.values, nu=5, C=1, debug=FALSE, filtering.iterations=50, print.iterations=100) {
   # old argments: permutation="ordered", which.lambda=1, 
   
   # parameters    
   #----------------------------------------------------------------------------- 
   aux   = starting.values    
   N     = dim(aux$Y)[2]     
   T     = dim(aux$Y)[1]     
   TT    = dim(aux$U)[1]     
   p     = T - TT
   M     = dim(aux$PR_TR)[1]
   d     = sum(restrictions$dj) 
   r     = ncol(restrictions$Q)
   
   # List of vectorized posteriors for each iteration    
   #-----------------------------------------------------------------------------     
   posteriors = list(        
      Sigma    = array(NA, c(N,N,M,S)),
      lambda.omega  = array(NA, c(N,M,S)),
      A        = array(NA, c(N,N,S)),
      a        = array(NA, c(r,S)),
      alpha    = array(NA, c(1+p*N,N,S)),
      PR_TR    = array(NA, c(M,M,S)),
      hyper    = array(NA, c(3,S)),
      w        = array(NA, c(d,S)),
      S.t      = array(NA, c(TT,S)),
      U.SF     = array(NA, c(dim(aux$U.SF),S)),
      MH.logkernel = rep(NA,S)
   )
   
   # S iterations of the Gibbs    
   #----------------------------------------------------------------------------- 
   for(iteration in 1:S){
      # Filtering and smoothing step      
      if(debug) print("Filtering and Smoothing step")
      tmp      = try(B.SVAR.MSH.identified.filtering.smoothing(aux, iterations=filtering.iterations))
      if(!inherits(tmp, "try-error")){aux = tmp}
      
      # Hidden Markov Chain step
      if(debug) print("Hidden Markov Chain step") 
      tmp      = try(B.SVAR.MSH.identified.hidden.markov.chain(aux, priors, restrictions))
      if(!inherits(tmp, "try-error")){aux = tmp}
      
      ### Hyper step
      if(debug) print("Hyper-parameters step")        
      tmp      = try(B.SVAR.MSH.identified.hyper(aux, priors))
      if(!inherits(tmp, "try-error")){aux = tmp}
      
      # Inverted-Gamma step
      if(debug) print("Inverted-Gamma step")        
      tmp      = try(B.SVAR.MSH.identified.inverted.gamma(aux, priors))
      if(!inherits(tmp, "try-error")){aux = tmp}
      
      # Structural step with MH algorithm              
      if(debug) print("Structural step")  
      tmp      = try(B.SVAR.MSH.identified.structural.MH(aux, priors, restrictions, nu=nu, C=C))
      if(!inherits(tmp, "try-error")){aux = tmp}
      
      # Regression step        
      if(debug) print("Regression step")
      tmp      = try(B.SVAR.MSH.identified.regression(aux, priors))
      if(!inherits(tmp, "try-error")){aux = tmp}
      
      # Save posteriors as vectors        
      for(regime in 1:M){
         posteriors$Sigma[,,regime,iteration] = aux$Sigma[,,regime]
      }
      posteriors$lambda.omega[,,iteration]   = as.matrix(aux$lambda.omega)
      posteriors$A[,,iteration]        = aux$A
      posteriors$a[,iteration]         = as.matrix(aux$a)
      posteriors$alpha[,,iteration]    = aux$alpha
      posteriors$PR_TR[,,iteration]    = aux$PR_TR
      posteriors$w[,iteration]         = matrix(aux$w, ncol=1)
      posteriors$hyper[,iteration]     = aux$hyper
      posteriors$S.t[,iteration]       = matrix(max.col(t(aux$xi)), ncol=1)
      posteriors$U.SF[,,iteration]        = aux$U.SF
      posteriors$MH.logkernel[iteration] = aux$MH.logkernel
      # Print iteration results
      if((debug) && ((iteration %% print.iterations)==0)){
         cat("---------------------------------------------------------- \nIteration:",iteration,"/", S,"\n")
         cat("Count transitions\n")
         print(count.regime.transitions(aux$xi))
         cat("lambda.omega\n")
         print(aux$lambda.omega)
         cat("A\n")
         print(aux$A)
      }else if((iteration %% print.iterations)==0) cat(" ",iteration)
   }
   
   cat("\n Structural step acceptance rate of MH algorithm: ", 1-rejectionRate(as.mcmc(posteriors$a[1,])))
   
   # Output    
   #-----------------------------------------------------------------------------    
   output  <- list(                
      last.draws  = aux,                
      posteriors  = posteriors               
   )
   
   return(output)
}
