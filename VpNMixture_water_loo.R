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

library(R2jags)
library(loo)
library(HDInterval)
library(MCMCvis)
library(ggplot2)
library(patchwork)
library(scales)
library(bayesplot)

dat <- read.csv("vaAllwater.csv", fill = FALSE, header = TRUE) 
vibrio <- list(c=dat$Vp.path,v=dat$Sample.Volume,samp=dat$FID,site=dat$Site.Num,temp=dat$temp.t,pheo=dat$pheo.t,turb=dat$turb.t,chlo=dat$chlo.t) #data string, total vibrio

# Base and Full N-Mixture model for serial dilution data 
cat(
  "model{
    for (i in 1:288) {
      # Observation model across serial dilutions
      c[i] ~ dbin(p[i],3)
      p[i] <- 1-exp(-MPN[i]*v[i])
      
      # likelihood function for loo
      log.like[i] <- log(pbin(c[i],p[i],3))

      # Biological model for microbial abundance
      MPN[i] ~ dpois(lambda[i])
      log(lambda[i]) <- b0 
      
    }
    b0 ~ dnorm(0,0.1)
  }",
  file="modelbase.jag"
)

cat(
  "model{
    for (i in 1:288) {
      # Observation model across serial dilutions
      c[i] ~ dbin(p[i],3)
      p[i] <- 1-exp(-MPN[i]*v[i])
      
      # Biological model for microbial abundance
      MPN[i] ~ dpois(lambda[i])
      log(lambda[i]) <- b0 + b1*temp[i] + b2*pheo[i] + b3*temp[i]*pheo[i] + U[site[i]] + V[samp[i]]
    }
    for (s in 1:3) {
      U[s] ~ dnorm(0,tau_U)
    }
    for (t in 1:48) {
      V[t] ~ dnorm(0,tau_V)
    }
    tau_U ~ dgamma(0.1,0.1)
    tau_V ~ dgamma(0.1,0.1)
    b0 ~ dnorm(0,0.1)
    b1 ~ dnorm(0,0.1)
    b2 ~ dnorm(0,0.1)
    b3 ~ dnorm(0,0.1)
  }",
  file="model1.jag"
)

# Initial params BOTH YEARS
base.inits <- list(list("b0"=1),
                   list("b0"=1),
                   list("b0"=1))

base.inits <- list(list("U"=numeric(3),"V"=numeric(48),"tau_U"=0.1,"tau_V"=0.1,"b0"=1,"b1"=0,"b2"=0,"b3"=0),
                   list("U"=numeric(3),"V"=numeric(48),"tau_U"=0.01,"tau_V"=0.1,"b0"=1,"b1"=0,"b2"=0,"b3"=0),
                   list("U"=numeric(3),"V"=numeric(48),"tau_U"=1,"tau_V"=0.1,"b0"=1,"b1"=0,"b2"=0,"b3"=0))

params <- c("b2")

m <- jags(data = vibrio,
               inits = base.inits,
               parameters.to.save = params,
               model.file = "model1.jag",
               n.chains = 3,
               n.iter = 10000,
               n.burnin = 1000,
               n.thin = 3)

m.parmlist <- m$BUGSoutput$sims.list
m.b2 <- m.parmlist$b2 
m.P <- 1 - length(m.b2[m.b2>0])/length(m.b2)
m

m.loo <- loo(m.loglike, r_eff = NA)
m.loo

m.mcmc <- as.mcmc(m)
m.gel <- gelman.diag(m.mcmc, confidence = 0.95, transform=FALSE, autoburnin=TRUE,
            multivariate=TRUE)
m.gel
# GR: Gelman Rubin convergence statistic
# ESS: Effective samples size of 1,000,000 MCMC iterations
# Loo 225.4 93.0
# Expected log predictive density -112.7 -46.5 
