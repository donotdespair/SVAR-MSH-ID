### In progress: This algo maybe will need to be extended 
### to order over states a chosen (which.lambda) structural shocks variance

#---------------------------------------------------------------------------------------------------
# Gibbs: Inverted Gamma step for BVAR models as in Markun (2012)
#---------------------------------------------------------------------------------------------------

B.SVAR.MSH.identified.inverted.gamma = function(aux, priors){
    
   # Setup constants     
   #-----------------------------------------------------------------------------    
   N  = dim(aux$Y)[2]     
   T  = dim(aux$Y)[1]     
   TT = dim(aux$U)[1]     
   p  = T - TT 

   # Draw lambda
   lambda.a       = priors$lambda.a + 2*TT
   
   lo.tmp         = aux$lambda.omega
   lo.tmp[,1]     = 1
   l.tmp          = rep(0,N)
   for (n in  1:N){
      for (m in 1:M){
         l.tmp[n] = l.tmp[n] + as.numeric(crossprod(aux$U.SF[which(aux$xi[m,]==1),n])/lo.tmp[n,m])
      }
   }
   lambda.b       = l.tmp + priors$lambda.b

   draw.chisq.l   = rchisq(N,lambda.a)
   draw.lambda    = lambda.b/draw.chisq.l
   
   # Draw omega
   draw        = matrix(0,N,M)
   for (n in  1:N){
      for (m in 2:M){
         T.m         = sum(aux$xi[m,])
         omega.a     = priors$omega.a + T.m
         
         l.tmp       = as.numeric(crossprod(aux$U.SF[which(aux$xi[m,]==1),n])/draw.lambda[n])
         omega.b     = l.tmp + priors$omega.b
         draw.chisq.o= rchisq(1,omega.a)
         draw[n,m]   = omega.b/draw.chisq.o
      }
   }
   
   draw[,1]          = draw.lambda

   # Store
   aux$lambda.omega  = draw
   for (regime in 1:M){
      if (regime ==1){
         lambda.m    = aux$lambda.omega[,regime]
      } else {
         lambda.m    = aux$lambda.omega[,regime] * aux$lambda.omega[,1]
      }
      aux$Sigma[,,regime]     = solve(aux$A) %*% diag(lambda.m) %*% solve(t(aux$A))
   }
    
   # Output 
   return(aux)
}
