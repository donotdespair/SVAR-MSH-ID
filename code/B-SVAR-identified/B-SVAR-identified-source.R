#---------------------------------------------------------------------------------------------------
# Gibbs algorithms
#---------------------------------------------------------------------------------------------------

path.Restricted.BVAR    <- "./code/B-SVAR-identified/"
# Paths
#-----------------------------------------------------------------------------------------------
path.GIBBS  <- paste(path.Restricted.BVAR, "Functions/", sep="")
# path.GIBBS  <- "./Functions/"

# Libraries
#-----------------------------------------------------------------------------------------------
#install.packages("mvtnorm", dependencies=TRUE)
#install.packages("MCMCpack", dependencies=TRUE)
#install.packages("tmvtnorm", dependencies=TRUE)

library("mvtnorm")
library("MCMCpack")
library("tmvtnorm") 
library("vars")
library("coda")                                                         
library("permute")
library("nor1mix")

source(paste(path.GIBBS , "00 B-SVAR-identified.R", sep=""))
source(paste(path.GIBBS , "00 B-SVAR-identified-General functions.R", sep=""))   
source(paste(path.GIBBS , "00 B-SVAR-identified-priors.R", sep=""))
source(paste(path.GIBBS , "01 B-SVAR-identified-Initialization.R", sep=""))   
source(paste(path.GIBBS , "04 B-SVAR-identified-Hyper.R", sep=""))
source(paste(path.GIBBS , "04 B-SVAR-identified-Inverted-Gamma.R", sep=""))
source(paste(path.GIBBS , "05 B-SVAR-identified-Structural-MH.R", sep=""))
source(paste(path.GIBBS , "06 B-SVAR-identified-Regression.R", sep=""))
source(paste(path.GIBBS , "08 B-SVAR-identified-MDD.R", sep=""))
source(paste(path.GIBBS , "09 B-SVAR-identified-SDDR.R", sep=""))
