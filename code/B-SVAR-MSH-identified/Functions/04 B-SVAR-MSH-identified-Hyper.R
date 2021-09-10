### In progress: This algo maybe will need to be extended 
### to order over states a chosen (which.lambda) structural shocks variance

#---------------------------------------------------------------------------------------------------
# Gibbs: Inverted Gamma step for BVAR models as in Markun (2012)
#---------------------------------------------------------------------------------------------------

B.SVAR.MSH.identified.hyper = function(aux, priors){
    
   # Setup constants     
   #-----------------------------------------------------------------------------    
   N  = dim(aux$Y)[2]     
   T  = dim(aux$Y)[1]     
   TT = dim(aux$U.SF)[1]     
   p  = T - TT 
   r  = ncol(restrictions$Q)
   
   # gamma a
   ga.a           = priors$hyper.a + r
   ga.b           = priors$hyper.b + as.numeric(crossprod(aux$a))
   draw.chisq.a   = rchisq(1,ga.a)
   draw.a         = ga.b/draw.chisq.a
   
   # gamma mu
   gm.a           = priors$hyper.a + N
   gm.b           = priors$hyper.b + as.numeric(crossprod(aux$alpha[1,]))
   draw.chisq.m   = rchisq(1,gm.a)
   draw.m         = gm.b/draw.chisq.m
   
   # gamma beta
   gb.a           = priors$hyper.a + p*(N^2)
   gb.tmp         = 0
   for (n in 1:N){
      gb.tmp      = gb.tmp + as.numeric(aux$alpha[-1,n] %*% priors$beta.H %*% aux$alpha[-1,n])
   }
   gb.b           = priors$hyper.b + gb.tmp
   draw.chisq.b   = rchisq(1,gb.a)
   draw.b         = gb.b/draw.chisq.b
   
   # Store
   aux$hyper      = as.matrix(c(draw.a,draw.m,draw.b))
    
   # Output 
   return(aux)
}
