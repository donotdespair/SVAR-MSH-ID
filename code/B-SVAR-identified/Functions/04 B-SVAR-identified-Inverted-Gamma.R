
B.SVAR.identified.inverted.gamma = function(aux, priors){
    
   # Setup constants     
   #-----------------------------------------------------------------------------    
   N  = dim(aux$Y)[2]     
   T  = dim(aux$Y)[1]     
   TT = dim(aux$U)[1]     
   p  = T - TT 

   # Draw lambda
   T.m            = TT
   lambda.a       = priors$lambda.a + 2*T.m
   
   l.tmp          = rep(0,N)
   for (n in  1:N){
      l.tmp[n] = l.tmp[n] + as.numeric(crossprod(aux$U.SF[,n]))
   }
   lambda.b       = l.tmp + priors$lambda.b

   draw.chisq.l   = rchisq(N,lambda.a)
   draw.lambda    = lambda.b/draw.chisq.l
   draw           = matrix(draw.lambda,ncol=1)

   # Store
   aux$lambda.omega  = draw
   aux$Sigma         = solve(aux$A) %*% diag(as.vector(draw)) %*% solve(t(aux$A))

   # Output 
   return(aux)
}
