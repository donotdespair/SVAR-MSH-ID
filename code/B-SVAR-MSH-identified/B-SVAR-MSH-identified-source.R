#---------------------------------------------------------------------------------------------------
# Gibbs algorithms
#---------------------------------------------------------------------------------------------------

path.Restricted.BVAR    <- "./code/B-SVAR-MSH-identified/"
# Paths
#-----------------------------------------------------------------------------------------------
path.GIBBS  <- paste(path.Restricted.BVAR, "Functions/", sep="")

# Libraries
#-----------------------------------------------------------------------------------------------
library("mvtnorm")
library("MCMCpack")
library("tmvtnorm") 
library("vars")
library("coda")                                                         
library("permute")
library("nor1mix")

source(paste(path.GIBBS , "00 B-SVAR-MSH-identified.R", sep=""))
source(paste(path.GIBBS , "00 B-SVAR-MSH-identified-General functions.R", sep=""))   
source(paste(path.GIBBS , "00 B-SVAR-MSH-identified-priors.R", sep=""))
source(paste(path.GIBBS , "01 B-SVAR-MSH-identified-Initialization.R", sep=""))   
source(paste(path.GIBBS , "02 B-SVAR-MSH-identified-Filtering-Smoothing.R", sep=""))   
source(paste(path.GIBBS , "03 B-SVAR-MSH-identified-Hidden-Markov-Chain.R", sep=""))   
source(paste(path.GIBBS , "04 B-SVAR-MSH-identified-Hyper.R", sep=""))
source(paste(path.GIBBS , "04 B-SVAR-MSH-identified-Inverted-Gamma.R", sep=""))
source(paste(path.GIBBS , "05 B-SVAR-MSH-identified-Structural-MH.R", sep=""))
source(paste(path.GIBBS , "06 B-SVAR-MSH-identified-Regression.R", sep=""))
source(paste(path.GIBBS , "08 B-SVAR-MSH-identified-MDD.R", sep=""))
source(paste(path.GIBBS , "09 B-SVAR-MSH-identified-SDDR.R", sep=""))
