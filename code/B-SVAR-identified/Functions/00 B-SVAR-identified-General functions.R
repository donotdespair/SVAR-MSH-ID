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


B.SVAR.identified.short.mcmc = function(qqq,keep.every=10){
   
   if (!is.null(qqq$posteriors$lambda.omega))      S = dim(qqq$posteriors$lambda.omega)[3]
   if (!is.null(qqq$posteriors$A))           S = dim(qqq$posteriors$A)[3]
   if (!is.null(qqq$posteriors$alpha))       S = dim(qqq$posteriors$alpha)[3]
   
   N  = dim(qqq$last.draws$Y)[2]
   qqq.new  = qqq
   ind= seq(from=1, to=S, by=keep.every)
   if (!is.null(qqq$posteriors$Sigma))       qqq.new$posteriors$Sigma      = qqq$posteriors$Sigma[,,ind]
   if (!is.null(qqq$posteriors$lambda.omega))      qqq.new$posteriors$lambda.omega     = qqq$posteriors$lambda.omega[,ind]
   if (!is.null(qqq$posteriors$A))           qqq.new$posteriors$A          = qqq$posteriors$A[,,ind]
   if (!is.null(qqq$posteriors$a))           qqq.new$posteriors$a          = qqq$posteriors$a[,ind]
   if (!is.null(qqq$posteriors$alpha))       qqq.new$posteriors$alpha      = qqq$posteriors$alpha[,,ind]
   if (!is.null(qqq$posteriors$hyper))           qqq.new$posteriors$hyper          = qqq$posteriors$hyper[,ind]
   if (!is.null(qqq$posteriors$U.SF))       qqq.new$posteriors$U.SF      = qqq$posteriors$U.SF[,,ind]
   if (!is.null(qqq$posteriors$MH.logkernel))       qqq.new$posteriors$MH.logkernel      = qqq$posteriors$MH.logkernel[ind]
   
   return(qqq.new)
}

B.SVAR.identified.merge = function(qqq1,qqq2){
   
   N  = dim(qqq1$last.draws$Y)[2]
   
   qqq.new     = qqq1
   qqq.new$posteriors$Sigma      = abind(qqq1$posteriors$Sigma,qqq2$posteriors$Sigma)
   qqq.new$posteriors$lambda.omega     = abind(qqq1$posteriors$lambda.omega,qqq2$posteriors$lambda.omega)
   qqq.new$posteriors$A          = abind(qqq1$posteriors$A,qqq2$posteriors$A)
   qqq.new$posteriors$a          = abind(qqq1$posteriors$a,qqq2$posteriors$a)
   qqq.new$posteriors$alpha      = abind(qqq1$posteriors$alpha,qqq2$posteriors$alpha)
   qqq.new$posteriors$hyper      = abind(qqq1$posteriors$hyper,qqq2$posteriors$hyper)
   qqq.new$posteriors$U.SF        = abind(qqq1$posteriors$U.SF,qqq2$posteriors$U.SF)
   qqq.new$posteriors$MH.logkernel        = c(qqq1$posteriors$MH.logkernel,qqq2$posteriors$MH.logkernel)
   return(qqq.new)
}
