#---------------------------------------------------------------------------------------------------
# EM algorithms
#---------------------------------------------------------------------------------------------------

file.path = "./code/EM-MSVAR/"

# Libraries
#---------------------------------------------------------------------------------------------------
library("mvtnorm")



# Functions common to all MSVAR models
#---------------------------------------------------------------------------------------------------
source(paste(file.path,"Functions/MSVAR-Expectation.R",sep="")) 

# Model specific functions
#---------------------------------------------------------------------------------------------------

# MSH 
source(paste(file.path,"Functions/MSH-Initialization.R",sep="")) 
source(paste(file.path,"Functions/MSH-Maximization.R",sep="")) 
source(paste(file.path,"Functions/MSH-VAR.R",sep="")) 
