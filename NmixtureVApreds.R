rm(list=ls(all=TRUE))

# Download JAGS at https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Windows/
# Load these packages the first time your run R
#install.packages("RTools")
#install.packages("rjags")
#install.packages("coda")
#install.packages("R2jags")
#install.packages("hdi")
#install.packages("MCMCvis")
#install.packages("patchwork")
#install.packages("ggdistribute")

library(R2jags)
library(loo)
library(HDInterval)
library(MCMCvis)
library(ggplot2)
library(patchwork)
library(scales)
library(ggdistribute)

dat <- read.csv("vaAllwater.csv", fill = FALSE, header = TRUE) 

# Base and Full N-Mixture model for serial dilution data 
cat(
  "model{
    for (i in 1:444) {
      # Observation model across serial dilutions
      c[i] ~ dbin(p[i],3)
      p[i] <- 1-exp(-MPN[samp[i]]*v[i])
      
      # Biological model for microbial abundance
    }
    for (s in 1:71) {
      MPN[s] ~ dgamma(0.1,0.1)
    }
  }",
  file="pred.jag"
)

# Initial params BOTH YEARS
pred.inits <- list(list("MPN"=numeric(71)+10),
                   list("MPN"=numeric(71)+5),
                   list("MPN"=numeric(71)+20))

parameters <- c("MPN")

vibrio <- list(c=dat$vvha,v=dat$volume,samp=dat$fid)

m.base <- jags(data = vibrio,
               inits = pred.inits,
               parameters.to.save = parameters,
               model.file = "pred.jag",
               n.chains = 3,
               n.iter = 10000,
               n.burnin = 1000,
               n.thin = 3)

out<-MCMCpstr(m.base,
              params = parameters,
              func = median,
              type = 'summary')

