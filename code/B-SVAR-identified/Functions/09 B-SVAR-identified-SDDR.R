# rm(list=ls())

# setwd("~/Dropbox/Codes/B-SVAR-MSH-identified")
# source("B-SVAR-MSH-identified-source.R")
# source("/Users/szeridan/Dropbox/Codes/EM-MSVAR/Caspur-MSVAR-source-EM.R")
# load("/Users/szeridan/Dropbox/Honey Bees/02 data BI2015/dataBI2015.RData")
# 
# setwd("/Users/szeridan/Dropbox/Honey Bees/03 empirical")
# load("classic-p4-M2-15.RData")
# qqq      = B.SVAR.MSH.identified.short.mcmc(qqq6,keep.every=1000)
# rm(qqq5,qqq6)



dratioig2     = function(x,a1,a2,b1,b2,logarithm=FALSE){
   
   lcratioig2  = (a1/2)*log(b1) + (a2/2)*log(b2) - lbeta((a1/2),(a2/2))
   lkratioig2  = ((a2-2)/2)*log(x) - ((a1+a2)/2)*log(b1 + b2*x)
   out      = as.numeric(lcratioig2 + lkratioig2)
   
   if (!logarithm) out = exp(out)
   return(out)
}

B.SVAR.identified.SDDR    = function(qqq,priors){
   aux         = qqq$last.draws
   posteriors  = qqq$posteriors
   M           = dim(aux$PR_TR)[1]
   p           = (dim(aux$Y)[1] - dim(aux$U)[1])
   T           = dim(aux$Y)[1]
   TT          = T - p
   N           = dim(aux$Y)[2] 
   S           = dim(posteriors$alpha)[3]
   
   log.prior      = matrix(NA, N, N)
   for (i in 1:N){
      for (j in 1:N){
         if (i!=j){
            log.prior[i,j]       = M*dratioig2(x=1, a1=priors$omega.a, a2=priors$omega.a, b1=priors$omega.b, b2=priors$omega.b, logarithm=TRUE)
         }
      }
   }
   
   posterior.omega.a    = priors$omega.a + apply(qqq$posteriors$S.t,2,count.regime.occurences.tmp)
   posterior.omega.b    = array(NA, c(N,M-1,S))
   
   log.posterior        = matrix(NA, N, N)
   log.posterior.tmp    = array(0, c(N, N, M-1, S))
   
   for (s in 1:S){
      for (m in 2:M){
         m.indicator          = qqq$posteriors$S.t[,s]==m
         posterior.omega.b[,m-1,s]    = priors$omega.b + (sum(m.indicator)-1)*apply(qqq$posteriors$U.SF[m.indicator,,s],2,var)/qqq$posteriors$lambda.omega[,1,s]
         for (i in 1:N){
            for (j in 1:N){
               if (i!=j){
                  log.posterior.tmp[i,j,m-1,s]   = dratioig2(
                     x = 1, 
                     a1 = posterior.omega.a[m,s], 
                     a2 = posterior.omega.a[m,s], 
                     b1 = posterior.omega.b[i,m-1,s], 
                     b2 = posterior.omega.b[j,m-1,s], 
                     logarithm = TRUE
                     )
               }
            }
         }
      }
   }
   
   log.posterior.s   = apply(log.posterior.tmp,c(1,2,4),sum)
   
   for (i in 1:N){
      for (j in 1:N){
         if (i!=j){
            lps                  = log.posterior.s[i,j,]
            mlps                 = max(lps)
            log.posterior[i,j]   = mlps + log(sum(exp(lps - mlps))) - log(length(lps))
    
         }
      }
   }
   
   log.SDDR    = log.posterior - log.prior
   
   output      = abind(log.SDDR,log.posterior,log.prior, along=3, new.names = c("log.SDDR","log.posterior","log.prior"))
   
   return(output)
}

count.regime.occurences.tmp = function(S.t){
   M           = max(S.t)
   occurences  = rep(NA,M)
   for (m in 1:M){
      occurences[m]  = sum(S.t==m)
   }
   return(occurences)
}
