#---------------------------------------------------------------------------------------------------
# Utility funstions for structural analysis SVECM:
# List of functions:
#     B.IRFfromSVECM(posterior,h)
#     B.VECM2VAR(posterior,priors)
#     VAR2VMA(q,vecA,N,p)
#     VMA2IRF(posterior.vma)
#     IRF2CumulativeIRF(posterior.irf)
#     IRF.lines(posterior.cirf, probability=.68)
#---------------------------------------------------------------------------------------------------





B.IRFfromSVECM = function(qqq,h) {
   # Requires normalized output
   posterior.var  = B.VAR(qqq)
   posterior.vma  = VAR2VMA(posterior.var,h)
   irf            = VMA2IRF(posterior.vma)
   cirf           = IRF2CumulativeIRF(irf)
   output         = list(
      irf            = irf,
      cumulative.irf = cirf
      )
   return(output)
}

B.IRFfromSVECM.standardized = function(qqq,h) {
  # Requires normalized output
  posterior.var  = B.VAR.standardized(qqq)
  posterior.vma  = VAR2VMA(posterior.var,h)
  irf            = VMA2IRF(posterior.vma)
  cirf           = IRF2CumulativeIRF(irf)
  output         = list(
    irf            = irf,
    cumulative.irf = cirf
  )
  return(output)
}
#---------------------------------------------------------------------------------------------------
B.VAR = function(posterior){
   aux         = posterior$last.draws
   posteriors  = posterior$posteriors
   
   # Setup constants     
   S  = dim(posteriors$Gamma)[3]
   N  = dim(aux$Y)[2]     
   T  = dim(aux$Y)[1]     
   TT = dim(aux$U)[1]     
   p  = T - TT
   
   # form output
   VAR.out        = new.env()
   VAR.out$A      = posteriors$A
   VAR.out$Beta   = posteriors$Gamma
   VAR.out        = as.list(VAR.out)
   for (s in 1:S){
      VAR.out$Beta[,,s]     = t(solve(posteriors$A[,,s])%*%t(VAR.out$Beta[,,s]))
   }
   # Output     
   return(VAR.out)
}

#---------------------------------------------------------------------------------------------------
B.VAR.standardized = function(posterior){
  aux         = posterior$last.draws
  posteriors  = posterior$posteriors
  
  # Setup constants     
  S  = dim(posteriors$Gamma)[3]
  N  = dim(aux$Y)[2]     
  T  = dim(aux$Y)[1]     
  TT = dim(aux$U)[1]     
  p  = T - TT
  
  # form output
  VAR.out        = new.env()
  VAR.out$A      = posteriors$A
  VAR.out$Beta   = posteriors$Gamma
  VAR.out        = as.list(VAR.out)
  for (s in 1:S){
    VAR.out$Beta[,,s]     = t(solve(posteriors$A[,,s])%*%t(VAR.out$Beta[,,s]))
    VAR.out$A[,,s]        = diag(diag(VAR.out$A[,,s])^(-1)) %*% VAR.out$A[,,s]
  }
  # Output     
  return(VAR.out)
}
#---------------------------------------------------------------------------------------------------
var4vma = function(q,Beta,N,p){
   # q      - order of the VMA process
   # vecA   - a N(1+pN)x1 vector of a vectorized 1+Np x N matrix A
   
   Beta     = t(Beta)
   B        = vector("list",p)
   for (i in 1:p){
      B[[i]]   = Beta[,2:(1+N)+(i-1)*N]
   }
   
   phi      = array(NA,c(N,N,q+1))
   phi[,,1] = diag(N)
   for (i in 1:q){
      PHI   = matrix(0,N,N)
      for (j in 1:min(p,i)){
         PHI      = PHI + phi[,,i+1-j]%*%B[[j]]
      }
      phi[,,i+1]  = PHI
   }
   return(phi)
}

#---------------------------------------------------------------------------------------------------
VAR2VMA = function(posterior.var,q){
   S        = dim(posterior.var$A)[3]
   N        = dim(posterior.var$A)[1]
   K        = dim(posterior.var$Beta)[1]
   p        = (K-1)/N
   Beta.post= posterior.var$Beta
   
   output   = array(NA,c(N,N,q+1,S))
   for (i in 1:S){
      output[,,,i]   = var4vma(q=q,Beta=Beta.post[,,i],N=N,p=p)
   }
   post.out = list(Phi = output, A = posterior.var$A)
   return(post.out)
}

VMA2IRF = function(posterior.vma){
   # Requires normalized output
  N     = dim(posterior.vma$Phi)[1]
  q     = dim(posterior.vma$Phi)[3]
  S     = dim(posterior.vma$Phi)[4]
  
  irf   = array(NA, dim(posterior.vma$Phi))
  for (i in 1:S){
    A.inv   = solve(posterior.vma$A[,,i])
    for (j in 1:q){
      irf[,,j,i] = posterior.vma$Phi[,,j,i]%*%A.inv
    }
  }
  return(irf)
}

IRF2CumulativeIRF = function(posterior.irf){
   # Requires normalized output
   N     = dim(posterior.irf)[1]
   h     = dim(posterior.irf)[3]
   S     = dim(posterior.irf)[4]
   cirf  = array(NA,dim(posterior.irf))
   for (i in 1:N){
      for (j in 1:N){
         for (s in 1:S){
            cirf[i,j,,s]      = cumsum(posterior.irf[i,j,,s])
         }
      }
   }
   return(cirf)
}

IRF.lines   = function(posterior.cirf, probability=.68){
   N     = dim(posterior.cirf$irf)[1]
   h     = dim(posterior.cirf$irf)[3]
   S     = dim(posterior.cirf$irf)[4]
   
   prob     = (1-probability)/2
   
   irf      = new.env()
   irf$upper = array(NA,dim(posterior.cirf$irf)[1:3])
   irf$irf   = array(NA,dim(posterior.cirf$irf)[1:3])
   irf$lower = array(NA,dim(posterior.cirf$irf)[1:3])
   irf      = as.list(irf)
   
   cirf     = new.env()
   cirf$upper = array(NA,dim(posterior.cirf$irf)[1:3])
   cirf$cirf  = array(NA,dim(posterior.cirf$irf)[1:3])
   cirf$lower = array(NA,dim(posterior.cirf$irf)[1:3])
   cirf     = as.list(cirf)
   
   for (i in 1:N){
      for (j in 1:N){
         for (en in 1:h){
            irf$lower[i,j,en]    = quantile(posterior.cirf$irf[i,j,en,], probs=prob)
            irf$irf[i,j,en]      = quantile(posterior.cirf$irf[i,j,en,], probs=.5)
            irf$upper[i,j,en]    = quantile(posterior.cirf$irf[i,j,en,], probs=1-prob)
            
            cirf$lower[i,j,en]   = quantile(posterior.cirf$cumulative.irf[i,j,en,], probs=prob)
            cirf$cirf[i,j,en]    = quantile(posterior.cirf$cumulative.irf[i,j,en,], probs=.5)
            cirf$upper[i,j,en]   = quantile(posterior.cirf$cumulative.irf[i,j,en,], probs=1-prob)
         }
      }
   }
   output = list(
      irf   = irf,
      cirf  = cirf
      )
   return(output)
}





# posterior.cirf = B.IRFfromSVECM(posterior=qqn,h=10)
# 
# line = IRF.lines(posterior.cirf,probability=.68)
# 
# line$irf$irf
# line$irf$lower
# line$irf$upper
# line$irf
# plot.irf.all(irf.line=line, name="barley", lines.width = 3)


# plot.irf.all = function(irf.line,name="", lines.width = 3){
#    N     = dim(irf.line$irf$irf)[1]
#    h     = dim(irf.line$irf$irf)[3]
#    
#    for (i in 1:N){
#       for (j in 1:N){
#          x  = 1:dim(irf.line$irf$irf)[3]
#          pdf(file=paste("./grphs/irf-",name,"-to-i",i,"-j",j,".pdf",sep=""),height=7,width=10)
#          par(mar=c(5,4,4,5)+.1)
#          plot(x,irf.line$irf$irf[i,j,], ylim=range(irf.line$irf),type="l",lwd=lines.width, col="black", xlab="time horizon",ylab="IRF",bty="n")
#          lines(irf.line$irf$lower[i,j,], col="black", lty=3,lwd=lines.width)
#          lines(irf.line$irf$upper[i,j,], col="black", lty=3,lwd=lines.width)
#          
# #          par(new=TRUE)
# #          plot(x,irf.line$cirf$cirf[i,j,], ylim=range(irf.line$cirf),type="l",lwd=lines.width, col="gray60",bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
# #          lines(irf.line$cirf$lower[i,j,], col="gray60", lty=3,lwd=lines.width)
# #          lines(irf.line$cirf$upper[i,j,], col="gray60", lty=3,lwd=lines.width)
# #          axis(4)
# #          mtext("Cumulative IRF",side=4,line=3)
#          # legend("topleft",col=c("black","black","gray60","gray60"),lty=c(1,3,1,3), lwd=rep(lines.width,4),legend=c("IRF","IRF band","CIRF","CIRF band"),horiz=TRUE,box.lwd=0)
#          dev.off()
#       }
#    }
# }




plot.irf.all = function(which.shock=1, irf.line,name="", lines.width = 3){
   N     = dim(irf.line$irf$irf)[1]
   h     = dim(irf.line$irf$irf)[3]
   j     = which.shock
   
   for (i in 1:N){
#       for (j in 1:N){
         x  = 1:dim(irf.line$irf$irf)[3]
         pdf(file=paste("./irf-",name,"-i",i,"-j",j,".pdf",sep=""),height=7,width=10)
         par(mar=c(5,4,4,5)+.1)
         plot(x,irf.line$irf$irf[i,j,], ylim=range(irf.line$irf),type="l",lwd=lines.width, col="black", xlab="time horizon",ylab="IRF",bty="n")
         lines(irf.line$irf$lower[i,j,], col="black", lty=3,lwd=lines.width)
         lines(irf.line$irf$upper[i,j,], col="black", lty=3,lwd=lines.width)
         
         par(new=TRUE)
         plot(x,irf.line$cirf$cirf[i,j,], ylim=range(irf.line$cirf),type="l",lwd=lines.width, col="gray60",bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
         lines(irf.line$cirf$lower[i,j,], col="gray60", lty=3,lwd=lines.width)
         lines(irf.line$cirf$upper[i,j,], col="gray60", lty=3,lwd=lines.width)
         axis(4)
         mtext("Cumulative IRF",side=4,line=3)
         # legend("topleft",col=c("black","black","gray60","gray60"),lty=c(1,3,1,3), lwd=rep(lines.width,4),legend=c("IRF","IRF band","CIRF","CIRF band"),horiz=TRUE,box.lwd=0)
         dev.off()
      }
#    }
}



one.plot.irf.all = function(irf.line,name="", lines.width = 3){
   N     = dim(irf.line$irf$irf)[1]
   h     = dim(irf.line$irf$irf)[3]
   
   pdf(file=paste("./irf-",name,".pdf",sep=""),height=7,width=10)
   for (i in 1:N){
      for (j in 1:N){
         x  = 1:dim(irf.line$irf$irf)[3]
         par(mar=c(5,4,4,5)+.1)
         plot(x,irf.line$irf$irf[i,j,], ylim=range(irf.line$irf),type="l",lwd=lines.width, col="black", xlab="time horizon",ylab="IRF",bty="n")
         lines(irf.line$irf$lower[i,j,], col="black", lty=3,lwd=lines.width)
         lines(irf.line$irf$upper[i,j,], col="black", lty=3,lwd=lines.width)
         
         par(new=TRUE)
         plot(x,irf.line$cirf$cirf[i,j,], ylim=range(irf.line$cirf),type="l",lwd=lines.width, col="gray60",bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
         lines(irf.line$cirf$lower[i,j,], col="gray60", lty=3,lwd=lines.width)
         lines(irf.line$cirf$upper[i,j,], col="gray60", lty=3,lwd=lines.width)
         axis(4)
         mtext("Cumulative IRF",side=4,line=3)
         # legend("topleft",col=c("black","black","gray60","gray60"),lty=c(1,3,1,3), lwd=rep(lines.width,4),legend=c("IRF","IRF band","CIRF","CIRF band"),horiz=TRUE,box.lwd=0)
      }
   }
   dev.off()
}




#---------------------------------------------------------------------------------------------------
# Structural.matrix = function(Sigma, restrictions){
#    # Restrict the structural matrix A according to the Rubio-Ramirez, Waggoner, Zha (2010) paper
#    
#    Sigma=matrix(qqq$posteriors$Sigma[,1000],ncol=4)
#    N     = dim(Sigma)[1]   
#    
#    restrictions = vector("list",N)
#    restrictions[[1]] = 0*diag(4)
#    restrictions[[2]] = rbind(diag(4)[1,], 0*diag(4)[2:4,])
#    restrictions[[3]] = rbind(diag(4)[1:2,], 0*diag(4)[3:4,])
#    restrictions[[4]] = rbind(diag(4)[1:3,], 0*diag(4)[4,])
#    A     = t(chol(solve(Sigma)))
#    A     = Sigma
#    P     = array(NA, dim(Sigma))
#    
#    for (i in 1:N){
#       if (i==1) {
#          P.i      = matrix(NA, 0, N)
#       } else {
#          P.i      = t(P[,1:(i-1)])
#       }      
#       Q.tilde.i   = rbind(restrictions[[i]] %*% A, P.i)
#       P[,i]       = qr.Q(qr(t(Q.tilde.i)))[,N]      
#    }
#    
#    final = A%*%P
#    
#    final[which(abs(final)<=1e-10)] = 0
#    final
#    
#    
#    
#    # 3-variate check
#    N     = 3
#    Sigma = rwishart(nu=5, V=diag(N))$W
#    A     = chol(S)
#    restrictions = vector("list",N)
#    restrictions[[1]]    = matrix(c(0,0,0,0,1,0,0,0,0),nrow=3)
#    restrictions[[2]]    = matrix(c(1,0,0,0,0,0,0,0,1),nrow=3)
#    restrictions[[3]]    = matrix(c(0,0,0,0,0,0,0,0,0),nrow=3)
#    solve(tcrossprod(A))
#    
#    
#    # M_i = Q_i %*% S.c
#    
#    # j=1
#    Q1.tilde = Q1%*%S.c
#    p1.t  = qr.Q(qr(t(Q1.tilde)))[,3]
#    
#    # j=2
#    Q2.tilde = rbind(Q2.bar%*%S.c,p1.t)
#    p2.t  = qr.Q(qr(t(Q2.tilde)))[,3]
#    
#    # j=3
#    Q3.tilde = rbind(Q3%*%S.c,p1.t,p2.t)
#    p3.t  = qr.Q(qr(t(Q3.tilde)))[,3]
#    
#    P = cbind(p1.t,p2.t,p3.t)
#    
#    final = S.c%*%P
#    
#    final[which(abs(final)<=0.1^10)] = 0
#    final
#    
#    
# }
# 
# is.zero = function(x){
#    output = prod(x==rep(0,length(x)))
#    return(output==1)
# }
# 
# no.zero.rows=function(x){
#    return(x[!apply(x,1,is.zero),])
# }
# 
# #---------------------------------------------------------------------------------------------------
# # VAR4SVAR = function(posterior.var, restrictions){
# #    
# # }
# 
# #---------------------------------------------------------------------------------------------------
# Switch.Beirut.with.Allepo.Barley = function(VAR.out){
#    switch.out = VAR.out
#    switch.A = function(A){
#       out = matrix(A,ncol=4)[c(1:3,5,4,6:7,9,8,10:11,13,12,14:15,17,16,18:19,21,20,22:23,25,24),c(1,2,4,3)]
#       return(out)
#    }
#    switch.out$A = apply(VAR.out$A,2,switch.A)
#    
#    switch.Sigma = function(Sigma){
#       out = matrix(Sigma,ncol=4)[c(1,2,4,3),c(1,2,4,3)]
#       return(out)
#    }
#    switch.out$Sigma = apply(VAR.out$Sigma,2,switch.Sigma)
#    return(switch.out)
# }
# 
# 
# 
# 
# 
