
B.SVAR.identified.Gibbs = function(S, priors, restrictions, starting.values, nu=5, C=1, debug=FALSE, print.iterations=100) {
   # parameters    
   #----------------------------------------------------------------------------- 
   aux   = starting.values    
   N     = dim(aux$Y)[2]     
   T     = dim(aux$Y)[1]     
   TT    = dim(aux$U)[1]     
   p     = T - TT
   r     = ncol(restrictions$Q)
   
   # List of vectorized posteriors for each iteration    
   #-----------------------------------------------------------------------------     
   posteriors = list(        
      Sigma    = array(NA, c(N,N,S)),
      lambda.omega  = array(NA, c(N,1,S)),
      A        = array(NA, c(N,N,S)),
      a        = array(NA, c(r,S)),
      alpha    = array(NA, c(1+p*N,N,S)),
      hyper    = array(NA, c(3,S)),
      U.SF     = array(NA, c(dim(aux$U.SF),S)),
      MH.logkernel = rep(NA,S)
   )
   
   # S iterations of the Gibbs    
   #----------------------------------------------------------------------------- 
   for(iteration in 1:S){
      ### Hyper step
      if(debug) print("Hyper-parameters step")        
      aux = B.SVAR.identified.hyper(aux, priors)  
      
      # Inverted-Gamma step
      if(debug) print("Inverted-Gamma step")        
      aux = B.SVAR.identified.inverted.gamma(aux, priors)
      
      # Structural step with MH algorithm              
      if(debug) print("Structural step")        
      aux = B.SVAR.identified.structural.MH(aux, priors, restrictions, nu=nu, C=C)
      
      # Regression step        
      if(debug) print("Regression step")        
      aux = B.SVAR.identified.regression(aux, priors)
      
      # Save posteriors as vectors        
      posteriors$Sigma[,,iteration]    = aux$Sigma
      posteriors$lambda.omega[,,iteration]   = as.matrix(aux$lambda.omega)
      posteriors$A[,,iteration]        = aux$A
      posteriors$a[,iteration]         = as.matrix(aux$a)
      posteriors$alpha[,,iteration]    = aux$alpha
      posteriors$hyper[,iteration]     = aux$hyper
      posteriors$U.SF[,,iteration]        = aux$U.SF
      posteriors$MH.logkernel[iteration] = aux$MH.logkernel
      # Print iteration results
      if((debug) && ((iteration %% print.iterations)==0)){
         cat("---------------------------------------------------------- \nIteration:",iteration,"/", S,"\n")
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
