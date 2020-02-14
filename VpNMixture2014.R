rm(list=ls(all=TRUE))

library(R2jags)
library(loo)
library(HDInterval)
library(MCMCvis)

dat <- read.csv("yr2seagrantvibrio.csv", fill = FALSE, header = TRUE) 
m1.dat<-list(c=dat$path,v=dat$mass,samp=dat$samp,gear=dat$gear,tide=dat$t2,mod=dat$mod,hi=dat$hi,time=dat$time)
