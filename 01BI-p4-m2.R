rm(list = ls()) 

options(echo=TRUE) 

source("./00codes/B-SVAR-MSH-identified/B-SVAR-MSH-identified-source.R")
  source("./00codes/EM-MSVAR/Caspur-MSVAR-source-EM.R")
load("dataBI2015.RData")

set.seed(1234567)

Y=100*Y
p  = 4
N  = dim(Y)[2]     
T  = dim(Y)[1]     
TT = T - p
M  = 2

# Belongia & Ireland (2015,JBES) Money in Taylor Rule
r  = 13
restrictions = list(
   Q  = matrix(0, N^2, r),
   q  = matrix(diag(N), ncol=1)
)

restrictions$Q[2:4,1:3]       = diag(3)
restrictions$Q[6,11]       = -1
restrictions$Q[9:11,4:6]     = diag(3)
restrictions$Q[21:21,7:7]   = diag(1)
restrictions$Q[24:24,8:8]   = diag(1)
restrictions$Q[27:28,9:10]   = diag(2)
restrictions$Q[30:30,11:11]   = diag(1)
restrictions$Q[33:33,12:12]   = diag(1)
restrictions$Q[35:35,13:13]   = diag(1)
restrictions$q[5,]            = -1

# restrictions$Q
# a     = rnorm(r)
# A0    = matrix(restrictions$Q%*%a + restrictions$q, ncol=N, nrow=N)

# For Transition probabilities
restrictions$M       = diag(M^2)
restrictions$dj      = seq(from=M, to=M, length.out=M)      

P                 = matrix(0, p*N, N)
P[1:N,]           = diag(N)
H.tmp             = rep(0,0)
for (i in 1:p) {H.tmp = c(H.tmp, rep(1/(i^2), N))}
H                 = diag(H.tmp)

priors   = list(
   beta.P         = P,
   beta.H         = H,
   lambda.a         = 1,
   lambda.b         = 1,
   omega.a         = 1,
   omega.b         = 3,
   hyper.a         = 1,
   hyper.b         = 1,
   w               = matrix(9 * diag(M) + matrix(1,M,M),ncol=1)    # Priors: Hidden Markov Chain, transition probabilities       
)

C = 0.1
# starting.values   = B.SVAR.MSH.VW.identified.initialization(Y, p, M, restrictions, which.Sigma=1)
# t.0      = proc.time()
# qqq      = B.SVAR.MSH.identified.Gibbs(S=1000, priors, restrictions, starting.values, debug=FALSE, print.iterations=100, C=C)
# t.1      = proc.time()  
# (t.1-t.0)/60
# 
# t.0      = proc.time()
# qqq1      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq$last.draws, debug=FALSE, print.iterations=100, C=C)
# t.1      = proc.time()  
# t.for.time = (t.1-t.0)/60
# 
# t.0      = proc.time()
# qqq2      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq1$last.draws, debug=FALSE, print.iterations=100, C=C)
# t.1      = proc.time()  
# t.for.time = (t.1-t.0)/60
# save(qqq1,qqq2,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/bi-p",p,"-M",M,"-01.RData", sep=""))
# 
# rm(qqq1)
# t.0      = proc.time()
# qqq3      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq2$last.draws, debug=FALSE, print.iterations=100, C=C)
# t.1      = proc.time()  
# t.for.time = (t.1-t.0)/60
# 
# rm(qqq2)
# t.0      = proc.time()
# qqq4      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq3$last.draws, debug=FALSE, print.iterations=100, C=C)
# t.1      = proc.time()  
# t.for.time = (t.1-t.0)/60
# save(qqq3,qqq4,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/bi-p",p,"-M",M,"-02.RData", sep=""))
# 
# rm(qqq3)
# t.0      = proc.time()
# qqq5      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq4$last.draws, debug=FALSE, print.iterations=100, C=C)
# t.1      = proc.time()  
# t.for.time = (t.1-t.0)/60
# 
# rm(qqq4)
# t.0      = proc.time()
# qqq6      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq5$last.draws, debug=FALSE, print.iterations=100, C=C)
# t.1      = proc.time()  
# t.for.time = (t.1-t.0)/60
# save(qqq5,qqq6,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/bi-p",p,"-M",M,"-03.RData", sep=""))


load(paste("./00estimation/bi-p",p,"-M",M,"-03.RData", sep=""))
rm(qqq5)
t.0      = proc.time()
qqq1      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq6$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60

rm(qqq6)
t.0      = proc.time()
qqq2      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq1$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60
save(qqq1,qqq2,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/bi-p",p,"-M",M,"-04.RData", sep=""))


rm(qqq1)
t.0      = proc.time()
qqq3      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq2$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60

rm(qqq2)
t.0      = proc.time()
qqq4      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq3$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60
save(qqq3,qqq4,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/bi-p",p,"-M",M,"-05.RData", sep=""))

rm(qqq3)
t.0      = proc.time()
qqq5      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq4$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60

rm(qqq4)
t.0      = proc.time()
qqq6      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq5$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60
save(qqq5,qqq6,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/bi-p",p,"-M",M,"-06.RData", sep=""))


rm(qqq5)
t.0      = proc.time()
qqq1      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq6$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60

rm(qqq6)
t.0      = proc.time()
qqq2      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq1$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60
save(qqq1,qqq2,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/bi-p",p,"-M",M,"-07.RData", sep=""))


rm(qqq1)
t.0      = proc.time()
qqq3      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq2$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60

rm(qqq2)
t.0      = proc.time()
qqq4      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq3$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60
save(qqq3,qqq4,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/bi-p",p,"-M",M,"-08.RData", sep=""))

rm(qqq3)
t.0      = proc.time()
qqq5      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq4$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60

rm(qqq4)
t.0      = proc.time()
qqq6      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq5$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60
save(qqq5,qqq6,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/bi-p",p,"-M",M,"-09.RData", sep=""))


rm(qqq5)
t.0      = proc.time()
qqq1      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq6$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60

rm(qqq6)
t.0      = proc.time()
qqq2      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq1$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60
save(qqq1,qqq2,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/bi-p",p,"-M",M,"-10.RData", sep=""))


rm(qqq1)
t.0      = proc.time()
qqq3      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq2$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60

rm(qqq2)
t.0      = proc.time()
qqq4      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq3$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60
save(qqq3,qqq4,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/bi-p",p,"-M",M,"-11.RData", sep=""))

rm(qqq3)
t.0      = proc.time()
qqq5      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq4$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60

rm(qqq4)
t.0      = proc.time()
qqq6      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq5$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60
save(qqq5,qqq6,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/bi-p",p,"-M",M,"-12.RData", sep=""))


rm(qqq5)
t.0      = proc.time()
qqq1      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq6$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60

rm(qqq6)
t.0      = proc.time()
qqq2      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq1$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60
save(qqq1,qqq2,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/bi-p",p,"-M",M,"-13.RData", sep=""))


rm(qqq1)
t.0      = proc.time()
qqq3      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq2$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60

rm(qqq2)
t.0      = proc.time()
qqq4      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq3$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60
save(qqq3,qqq4,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/bi-p",p,"-M",M,"-14.RData", sep=""))

rm(qqq3)
t.0      = proc.time()
qqq5      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq4$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60

rm(qqq4)
t.0      = proc.time()
qqq6      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq5$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60
save(qqq5,qqq6,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/bi-p",p,"-M",M,"-15.RData", sep=""))


# 
# rm(qqq5)
# t.0      = proc.time()
# qqq1      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq6$last.draws, debug=FALSE, print.iterations=100, C=C)
# t.1      = proc.time()  
# t.for.time = (t.1-t.0)/60
# 
# rm(qqq6)
# t.0      = proc.time()
# qqq2      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq1$last.draws, debug=FALSE, print.iterations=100, C=C)
# t.1      = proc.time()  
# t.for.time = (t.1-t.0)/60
# save(qqq1,qqq2,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/bi-p",p,"-M",M,"-16.RData", sep=""))
# 
# 
# rm(qqq1)
# t.0      = proc.time()
# qqq3      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq2$last.draws, debug=FALSE, print.iterations=100, C=C)
# t.1      = proc.time()  
# t.for.time = (t.1-t.0)/60
# 
# rm(qqq2)
# t.0      = proc.time()
# qqq4      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq3$last.draws, debug=FALSE, print.iterations=100, C=C)
# t.1      = proc.time()  
# t.for.time = (t.1-t.0)/60
# save(qqq3,qqq4,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/bi-p",p,"-M",M,"-17.RData", sep=""))
# 
# rm(qqq3)
# t.0      = proc.time()
# qqq5      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq4$last.draws, debug=FALSE, print.iterations=100, C=C)
# t.1      = proc.time()  
# t.for.time = (t.1-t.0)/60
# 
# rm(qqq4)
# t.0      = proc.time()
# qqq6      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq5$last.draws, debug=FALSE, print.iterations=100, C=C)
# t.1      = proc.time()  
# t.for.time = (t.1-t.0)/60
# save(qqq5,qqq6,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/bi-p",p,"-M",M,"-18.RData", sep=""))
# 
# 
# 
# 
# rm(qqq5)
# t.0      = proc.time()
# qqq1      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq6$last.draws, debug=FALSE, print.iterations=100, C=C)
# t.1      = proc.time()  
# t.for.time = (t.1-t.0)/60
# 
# rm(qqq6)
# t.0      = proc.time()
# qqq2      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq1$last.draws, debug=FALSE, print.iterations=100, C=C)
# t.1      = proc.time()  
# t.for.time = (t.1-t.0)/60
# save(qqq1,qqq2,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/bi-p",p,"-M",M,"-19.RData", sep=""))
# 
# 
# rm(qqq1)
# t.0      = proc.time()
# qqq3      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq2$last.draws, debug=FALSE, print.iterations=100, C=C)
# t.1      = proc.time()  
# t.for.time = (t.1-t.0)/60
# 
# rm(qqq2)
# t.0      = proc.time()
# qqq4      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq3$last.draws, debug=FALSE, print.iterations=100, C=C)
# t.1      = proc.time()  
# t.for.time = (t.1-t.0)/60
# save(qqq3,qqq4,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/bi-p",p,"-M",M,"-20.RData", sep=""))
# 
# rm(qqq3)
# t.0      = proc.time()
# qqq5      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq4$last.draws, debug=FALSE, print.iterations=100, C=C)
# t.1      = proc.time()  
# t.for.time = (t.1-t.0)/60
# 
# rm(qqq4)
# t.0      = proc.time()
# qqq6      = B.SVAR.MSH.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq5$last.draws, debug=FALSE, print.iterations=100, C=C)
# t.1      = proc.time()  
# t.for.time = (t.1-t.0)/60
# save(qqq5,qqq6,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/bi-p",p,"-M",M,"-21.RData", sep=""))



























# p=4;M=2
# model = "bi"
# load(paste(model,"-p",p,"-M",M,"-13.RData", sep=""))
# qqq      = B.SVAR.MSH.identified.merge(B.SVAR.MSH.identified.short.mcmc(qqq1,keep.every=10),B.SVAR.MSH.identified.short.mcmc(qqq2,keep.every=10))
# rm(qqq1,qqq2)
# load(paste(model,"-p",p,"-M",M,"-14.RData", sep=""))
# qqq      = B.SVAR.MSH.identified.merge(qqq,B.SVAR.MSH.identified.short.mcmc(qqq3,keep.every=10))
# qqq      = B.SVAR.MSH.identified.merge(qqq,B.SVAR.MSH.identified.short.mcmc(qqq4,keep.every=10))
# rm(qqq3,qqq4)
# load(paste(model,"-p",p,"-M",M,"-15.RData", sep=""))
# qqq      = B.SVAR.MSH.identified.merge(qqq,B.SVAR.MSH.identified.short.mcmc(qqq5,keep.every=10))
# qqq      = B.SVAR.MSH.identified.merge(qqq,B.SVAR.MSH.identified.short.mcmc(qqq6,keep.every=10))
# rm(qqq5,qqq6)
# 
# 
# qqq      = B.SVAR.MSH.identified.short.mcmc(qqq,keep.every=1)
# 
# 
# t0       = proc.time()
# log.mdd  = B.SVAR.MSH.identified.MDD(Sp=100000,Gibbs.output=qqq, priors, restrictions)
# t1       = proc.time()
# t.for.time = (t1-t0)[3]
# save(log.mdd,qqq,priors,restrictions,t.for.time, file = paste("./00estimation/mdd-",model,"-p",p,"-M",M,".RData", sep=""))

# load(paste("./00estimation/mdd-",model,"-p",p,"-M",M,".RData", sep=""))
# log.sddr    = B.SVAR.MSH.identified.SDDR(qqq=qqq,priors)
# save(log.sddr,qqq,priors, file = paste("./00estimation/sddr-",model,"-p",p,"-M",M,".RData", sep=""))

