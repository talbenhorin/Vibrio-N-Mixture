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

dat <- read.csv("vaoysterwater.csv", fill = FALSE, header = TRUE) 
vibrio <- list(c=dat$path,v=dat$mass,samp=dat$fid,site=dat$site,temp=dat$stan.temp,water=dat$water,watervib=dat$stan.path) #data string, total vibrio

cat(
  "model{
    for (i in 1:288) {
      # Observation model across serial dilutions
      c[i] ~ dbin(p[i],3)
      p[i] <- 1-exp(-MPN[i]*v[i])
      
      # Biological model for microbial abundance
      MPN[i] ~ dpois(lambda[i])
      log(lambda[i]) <- b0 
      
      # likelihood function for loo
      like[i] <- log(pbin(c[i],p[i],3))
      
    }
    b0 ~ dnorm(0,0.1)

    #sum log likihoods for loo
    log.like <- sum(like)
    
  }",
  file="base.jag"
)

cat(
  "model{
    for (i in 1:288) {
      # Observation model across serial dilutions
      c[i] ~ dbin(p[i],3)
      p[i] <- 1-exp(-MPN[i]*v[i])
      
      # Biological model for microbial abundance
      MPN[i] ~ dpois(lambda[i])
      log(lambda[i]) <- b0 + U[site[i]] + V[samp[i]]
      
      # likelihood function for loo
      like[i] <- log(pbin(c[i],p[i],3))
      
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

    #sum log likihoods for loo
    log.like <- sum(like)
    
  }",
  file="m0.jag"
)

cat(
  "model{
    for (i in 1:288) {
      # Observation model across serial dilutions
      c[i] ~ dbin(p[i],3)
      p[i] <- 1-exp(-MPN[i]*v[i])
      
      # Biological model for microbial abundance
      MPN[i] ~ dpois(lambda[i])
      log(lambda[i]) <- b0 + b1*watervib[water[i]] + U[site[i]] + V[samp[i]]
      
      # likelihood function for loo
      like[i] <- log(pbin(c[i],p[i],3))
      
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

    #sum log likihoods for loo
    log.like <- sum(like)
    
  }",
  file="m1.jag"
)

cat(
  "model{
    for (i in 1:288) {
      # Observation model across serial dilutions
      c[i] ~ dbin(p[i],3)
      p[i] <- 1-exp(-MPN[i]*v[i])
      
      # Biological model for microbial abundance
      MPN[i] ~ dpois(lambda[i])
      log(lambda[i]) <- b0 + b1*watervib[water[i]] + b2*temp[i] + U[site[i]] + V[samp[i]]
      
      # likelihood function for loo
      like[i] <- log(pbin(c[i],p[i],3))
      
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

    #sum log likihoods for loo
    log.like <- sum(like)
    
  }",
  file="m2.jag"
)

cat(
  "model{
    for (i in 1:288) {
      # Observation model across serial dilutions
      c[i] ~ dbin(p[i],3)
      p[i] <- 1-exp(-MPN[i]*v[i])
      
      # Biological model for microbial abundance
      MPN[i] ~ dpois(lambda[i])
      log(lambda[i]) <- b0 + b1*watervib[water[i]] + b2*temp[i] + b3*watervib[water[i]]*watervib[water[i]] + U[site[i]] + V[samp[i]]
      
      # likelihood function for loo
      like[i] <- log(pbin(c[i],p[i],3))
      
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

    #sum log likihoods for loo
    log.like <- sum(like)
    
  }",
  file="m3.jag"
)

inits.base <- list(list("b0"=1),
                   list("b0"=0.5),
                   list("b0"=2))
inits.m0 <- list(list("U"=numeric(3),"V"=numeric(48),"tau_U"=0.1,"tau_V"=0.1,"b0"=1),
                   list("U"=numeric(3),"V"=numeric(48),"tau_U"=0.01,"tau_V"=0.1,"b0"=1),
                   list("U"=numeric(3),"V"=numeric(48),"tau_U"=1,"tau_V"=0.1,"b0"=1))
inits.m1 <- list(list("U"=numeric(3),"V"=numeric(48),"tau_U"=0.1,"tau_V"=0.1,"b0"=1,"b1"=0),
                   list("U"=numeric(3),"V"=numeric(48),"tau_U"=0.01,"tau_V"=0.1,"b0"=1,"b1"=0),
                   list("U"=numeric(3),"V"=numeric(48),"tau_U"=1,"tau_V"=0.1,"b0"=1,"b1"=0))
inits.m2 <- list(list("U"=numeric(3),"V"=numeric(48),"tau_U"=0.1,"tau_V"=0.1,"b0"=1,"b1"=0,"b2"=0),
                   list("U"=numeric(3),"V"=numeric(48),"tau_U"=0.01,"tau_V"=0.1,"b0"=1,"b1"=0,"b2"=0),
                   list("U"=numeric(3),"V"=numeric(48),"tau_U"=1,"tau_V"=0.1,"b0"=1,"b1"=0,"b2"=0))
inits.m3 <- list(list("U"=numeric(3),"V"=numeric(48),"tau_U"=0.1,"tau_V"=0.1,"b0"=1,"b1"=0,"b2"=0,"b3"=0),
                 list("U"=numeric(3),"V"=numeric(48),"tau_U"=0.01,"tau_V"=0.1,"b0"=1,"b1"=0,"b2"=0,"b3"=0),
                 list("U"=numeric(3),"V"=numeric(48),"tau_U"=1,"tau_V"=0.1,"b0"=1,"b1"=0,"b2"=0,"b3"=0))

params.base <- c("log.like","b0")
params.m0 <- c("log.like","b0")
params.m1 <- c("log.like","b0","b1")
params.m2 <- c("log.like","b0","b1","b2")
params.m3 <- c("log.like","b0","b1","b2","b3")

m <- jags(data = vibrio,
          inits = inits.m3,
          parameters.to.save = params.m3,
          model.file = "m3.jag",
          n.chains = 3,
          n.iter = 50000,
          n.burnin = 5000,
          n.thin = 3)

m.parmlist <- m$BUGSoutput$sims.list
m.loglike <- m.parmlist$log.like

m.loo <- loo(m.loglike, r_eff = NA)
m.loo
m.mcmc <- as.mcmc(m)
m.gel <- gelman.diag(m.mcmc, confidence = 0.95, transform=FALSE, autoburnin=TRUE,
                     multivariate=TRUE)
m.gel