rm(list = ls()) 

options(echo=TRUE) 

# Estimation of homoskedastic model with a lower-triangular matrix A_0 being estimated

source("./00codes/B-SVAR-identified/B-SVAR-identified-source.R")
load("dataBI2015.RData")

set.seed(1234567)

Y=100*Y
p  = 4
N  = dim(Y)[2]     
T  = dim(Y)[1]     
TT = T - p
model = "classic"
M  = 1

# lower-triangular matrix with ones on the diagonal
r  = 15
restrictions = list(
   Q  = matrix(0, N^2, r),
   q  = matrix(diag(N), ncol=1)
)

restrictions$Q[2:6,1:5]       = diag(5)
restrictions$Q[9:12,6:9]     = diag(4)
restrictions$Q[16:18,10:12]   = diag(3)
restrictions$Q[23:24,13:14]   = diag(2)
restrictions$Q[30:30,15:15]   = 1

# restrictions$Q
# a     = rnorm(15)
# A0    = matrix(restrictions$Q%*%a + restrictions$q, ncol=N, nrow=N)

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
   hyper.a         = 1,
   hyper.b         = 1
)

starting.values   = B.SVAR.identified.initialization(Y, p, restrictions)
C = 0.2
t.0      = proc.time()
qqq      = B.SVAR.identified.Gibbs(S=1000, priors, restrictions, starting.values, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
(t.1-t.0)/60

t.0      = proc.time()
qqq1      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60

t.0      = proc.time()
qqq2      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq1$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60
save(qqq1,qqq2,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/classic-p",p,"-M",M,"-01.RData", sep=""))

rm(qqq1)
t.0      = proc.time()
qqq3      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq2$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60

rm(qqq2)
t.0      = proc.time()
qqq4      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq3$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60
save(qqq3,qqq4,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/classic-p",p,"-M",M,"-02.RData", sep=""))

rm(qqq3)
t.0      = proc.time()
qqq5      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq4$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60

rm(qqq4)
t.0      = proc.time()
qqq6      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq5$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60
save(qqq5,qqq6,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/classic-p",p,"-M",M,"-03.RData", sep=""))

rm(qqq5)
t.0      = proc.time()
qqq1      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq6$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60

rm(qqq6)
t.0      = proc.time()
qqq2      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq1$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60
save(qqq1,qqq2,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/classic-p",p,"-M",M,"-04.RData", sep=""))


rm(qqq1)
t.0      = proc.time()
qqq3      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq2$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60

rm(qqq2)
t.0      = proc.time()
qqq4      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq3$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60
save(qqq3,qqq4,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/classic-p",p,"-M",M,"-05.RData", sep=""))

rm(qqq3)
t.0      = proc.time()
qqq5      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq4$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60

rm(qqq4)
t.0      = proc.time()
qqq6      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq5$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60
save(qqq5,qqq6,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/classic-p",p,"-M",M,"-06.RData", sep=""))


rm(qqq5)
t.0      = proc.time()
qqq1      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq6$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60

rm(qqq6)
t.0      = proc.time()
qqq2      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq1$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60
save(qqq1,qqq2,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/classic-p",p,"-M",M,"-07.RData", sep=""))


rm(qqq1)
t.0      = proc.time()
qqq3      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq2$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60

rm(qqq2)
t.0      = proc.time()
qqq4      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq3$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60
save(qqq3,qqq4,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/classic-p",p,"-M",M,"-08.RData", sep=""))

rm(qqq3)
t.0      = proc.time()
qqq5      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq4$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60

rm(qqq4)
t.0      = proc.time()
qqq6      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq5$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60
save(qqq5,qqq6,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/classic-p",p,"-M",M,"-09.RData", sep=""))



rm(qqq5)
t.0      = proc.time()
qqq1      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq6$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60

rm(qqq6)
t.0      = proc.time()
qqq2      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq1$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60
save(qqq1,qqq2,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/classic-p",p,"-M",M,"-10.RData", sep=""))


rm(qqq1)
t.0      = proc.time()
qqq3      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq2$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60

rm(qqq2)
t.0      = proc.time()
qqq4      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq3$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60
save(qqq3,qqq4,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/classic-p",p,"-M",M,"-11.RData", sep=""))

rm(qqq3)
t.0      = proc.time()
qqq5      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq4$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60

rm(qqq4)
t.0      = proc.time()
qqq6      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq5$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60
save(qqq5,qqq6,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/classic-p",p,"-M",M,"-12.RData", sep=""))



rm(qqq5)
t.0      = proc.time()
qqq1      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq6$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60

rm(qqq6)
t.0      = proc.time()
qqq2      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq1$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60
save(qqq1,qqq2,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/classic-p",p,"-M",M,"-13.RData", sep=""))


rm(qqq1)
t.0      = proc.time()
qqq3      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq2$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60

rm(qqq2)
t.0      = proc.time()
qqq4      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq3$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60
save(qqq3,qqq4,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/classic-p",p,"-M",M,"-14.RData", sep=""))

rm(qqq3)
t.0      = proc.time()
qqq5      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq4$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60

rm(qqq4)
t.0      = proc.time()
qqq6      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq5$last.draws, debug=FALSE, print.iterations=100, C=C)
t.1      = proc.time()  
t.for.time = (t.1-t.0)/60
save(qqq5,qqq6,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/classic-p",p,"-M",M,"-15.RData", sep=""))
# 
# 
# 
# rm(qqq5)
# t.0      = proc.time()
# qqq1      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq6$last.draws, debug=FALSE, print.iterations=1000)
# t.1      = proc.time()  
# t.for.time = (t.1-t.0)/60
# 
# rm(qqq6)
# t.0      = proc.time()
# qqq2      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq1$last.draws, debug=FALSE, print.iterations=1000)
# t.1      = proc.time()  
# t.for.time = (t.1-t.0)/60
# save(qqq1,qqq2,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/",model,"-p",p,"-M",M,"-16.RData", sep=""))
# 
# 
# rm(qqq1)
# t.0      = proc.time()
# qqq3      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq2$last.draws, debug=FALSE, print.iterations=1000)
# t.1      = proc.time()  
# t.for.time = (t.1-t.0)/60
# 
# rm(qqq2)
# t.0      = proc.time()
# qqq4      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq3$last.draws, debug=FALSE, print.iterations=1000)
# t.1      = proc.time()  
# t.for.time = (t.1-t.0)/60
# save(qqq3,qqq4,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/",model,"-p",p,"-M",M,"-17.RData", sep=""))
# 
# rm(qqq3)
# t.0      = proc.time()
# qqq5      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq4$last.draws, debug=FALSE, print.iterations=1000)
# t.1      = proc.time()  
# t.for.time = (t.1-t.0)/60
# 
# rm(qqq4)
# t.0      = proc.time()
# qqq6      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq5$last.draws, debug=FALSE, print.iterations=1000)
# t.1      = proc.time()  
# t.for.time = (t.1-t.0)/60
# save(qqq5,qqq6,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/",model,"-p",p,"-M",M,"-18.RData", sep=""))
# 
# #
# rm(qqq5)
# t.0      = proc.time()
# qqq7      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq6$last.draws, debug=FALSE, print.iterations=1000)
# t.1      = proc.time()  
# t.for.time = (t.1-t.0)/60
# 
# rm(qqq6)
# t.0      = proc.time()
# qqq8      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq7$last.draws, debug=FALSE, print.iterations=1000)
# t.1      = proc.time()  
# t.for.time = (t.1-t.0)/60
# save(qqq7,qqq8,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/",model,"-p",p,"-M",M,"-19.RData", sep=""))
# 
# 
# rm(qqq7)
# t.0      = proc.time()
# qqq9      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq8$last.draws, debug=FALSE, print.iterations=1000)
# t.1      = proc.time()  
# t.for.time = (t.1-t.0)/60
# 
# rm(qqq8)
# t.0      = proc.time()
# qqq10      = B.SVAR.identified.Gibbs(S=200000, priors, restrictions, starting.values=qqq9$last.draws, debug=FALSE, print.iterations=1000)
# t.1      = proc.time()  
# t.for.time = (t.1-t.0)/60
# save(qqq9,qqq10,priors,restrictions,starting.values,Y,t.for.time, file = paste("./00estimation/",model,"-p",p,"-M",M,"-20.RData", sep=""))







load(paste("./00estimation/",model,"-p",p,"-M",M,"-13.RData", sep=""))
qqq      = B.SVAR.identified.merge(B.SVAR.identified.short.mcmc(qqq1,keep.every=10),B.SVAR.identified.short.mcmc(qqq2,keep.every=10))
rm(qqq1,qqq2)
load(paste("./00estimation/",model,"-p",p,"-M",M,"-14.RData", sep=""))
qqq      = B.SVAR.identified.merge(qqq,B.SVAR.identified.short.mcmc(qqq3,keep.every=10))
qqq      = B.SVAR.identified.merge(qqq,B.SVAR.identified.short.mcmc(qqq4,keep.every=10))
rm(qqq3,qqq4)
load(paste("./00estimation/",model,"-p",p,"-M",M,"-15.RData", sep=""))
qqq      = B.SVAR.identified.merge(qqq,B.SVAR.identified.short.mcmc(qqq5,keep.every=10))
qqq      = B.SVAR.identified.merge(qqq,B.SVAR.identified.short.mcmc(qqq6,keep.every=10))
rm(qqq5,qqq6)
cat("\n Computing the lnMDD\n")
t0       = proc.time()
log.mdd  = B.SVAR.identified.MDD(Sp=100000,Gibbs.output=qqq, priors, restrictions)
t1       = proc.time()
t.for.time = (t1-t0)[3]
save(log.mdd,qqq,priors,restrictions,t.for.time, file = paste("./00estimation/",model,"-p",p,"-M",M,".RData", sep=""))





































 
# p=4;M=2
# model = "classic"
# load(paste(model,"-p",p,"-M",M,"-13.RData", sep=""))
# qqq      = B.SVAR.identified.merge(B.SVAR.identified.short.mcmc(qqq1,keep.every=10),B.SVAR.identified.short.mcmc(qqq2,keep.every=10))
# rm(qqq1,qqq2)
# load(paste(model,"-p",p,"-M",M,"-14.RData", sep=""))
# qqq      = B.SVAR.identified.merge(qqq,B.SVAR.identified.short.mcmc(qqq3,keep.every=10))
# qqq      = B.SVAR.identified.merge(qqq,B.SVAR.identified.short.mcmc(qqq4,keep.every=10))
# rm(qqq3,qqq4)
# load(paste(model,"-p",p,"-M",M,"-15.RData", sep=""))
# qqq      = B.SVAR.identified.merge(qqq,B.SVAR.identified.short.mcmc(qqq5,keep.every=10))
# qqq      = B.SVAR.identified.merge(qqq,B.SVAR.identified.short.mcmc(qqq6,keep.every=10))
# rm(qqq5,qqq6)
# 
# 
# qqq      = B.SVAR.identified.short.mcmc(qqq,keep.every=1)
# 
# 
# 
# t0       = proc.time()
# log.mdd  = B.SVAR.identified.MDD(Sp=100000,Gibbs.output=qqq, priors, restrictions)
# t1       = proc.time()
# t.for.time = (t1-t0)[3]
# save(log.mdd,qqq,priors,restrictions,t.for.time, file = paste("./00estimation/mdd-",model,"-p",p,"-M",M,".RData", sep=""))

# load(paste("./00estimation/mdd-",model,"-p",p,"-M",M,".RData", sep=""))
# log.sddr    = B.SVAR.identified.SDDR(qqq=qqq,priors)
# save(log.sddr,qqq,priors, file = paste("./00estimation/sddr-",model,"-p",p,"-M",M,".RData", sep=""))



# head(Y)
# qqq = qqq1
# plot.ts(qqq$posteriors$MH.logkernel)
# 1-rejectionRate(as.mcmc(qqq$posteriors$a[1,]))
# 
# plot.ts(t(qqq$posteriors$a[1:10,]))
# plot.ts(t(qqq$posteriors$a[11:20,]))
# plot.ts(t(qqq$posteriors$a[21:30,]))
# 
# plot.ts(t(qqq$posteriors$lambda.omega[,1,]))
# plot.ts(t(qqq$posteriors$lambda.omega[,2,]))
# plot.ts(t(qqq$posteriors$lambda.omega[,3,]))
# apply(qqq$posteriors$lambda.omega,1:2,mean)
# 
# plot.ts(t(qqq$posteriors$PR_TR[1,,]))
# plot.ts(t(qqq$posteriors$PR_TR[2,,]))
# plot.ts(t(qqq$posteriors$PR_TR[3,,]))
# apply(qqq$posteriors$PR_TR,1:2,mean)
# 
# plot.ts(t(qqq$posteriors$hyper))
# apply(qqq$posteriors$hyper,1,mean)
# 
# plot.ts(t(qqq$posteriors$alpha[1,,]))
# apply(qqq$posteriors$alpha,1:2,mean)
# 
# plot.ts(G.probabilities.regime(posteriors=qqq$posteriors, keep.iterations=1:100000))
