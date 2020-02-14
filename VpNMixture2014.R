rm(list=ls(all=TRUE))

library(R2jags)
library(loo)
library(HDInterval)
library(MCMCvis)

dat <- read.csv("seagrantvibrio.csv", fill = FALSE, header = TRUE) 

