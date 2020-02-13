rm(list=ls(all=TRUE))

# Download JAGS at https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Windows/
# Load these packages the first time your run R
#install.packages("RTools")
#install.packages("rjags")
#install.packages("coda")
#install.packages("R2jags")
#install.packages("hdi")
#install.packages("MCMCvis")

library(R2jags)
library(loo)
library(HDInterval)
library(MCMCvis)

dat <- read.csv("vaoysters.csv", fill = FALSE, header = TRUE) 
