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

dat <- read.csv("vaoysterwaterpilf.csv", fill = FALSE, header = TRUE) 
vibrio <- list(c=dat$pilf,v=dat$mass,samp=dat$new.fid,site=dat$site,temp=dat$stan.temp,water=dat$new.water,watervib=dat$stan.vvha) #data string, total vibrio

cat(
  "model{
    for (i in 1:1362) {
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
    for (i in 1:1362) {
      # Observation model across serial dilutions
      c[i] ~ dbin(p[i],3)
      p[i] <- 1-exp(-MPN[i]*v[i])
      
      # Biological model for microbial abundance
      MPN[i] ~ dpois(lambda[i])
      log(lambda[i]) <- b0 + U[site[i]] + V[samp[i]]
      
      # likelihood function for loo
      like[i] <- log(pbin(c[i],p[i],3))
      
    }
    for (s in 1:4) {
      U[s] ~ dnorm(0,tau_U)
    }
    for (t in 1:222) {
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
    for (i in 1:1362) {
      # Observation model across serial dilutions
      c[i] ~ dbin(p[i],3)
      p[i] <- 1-exp(-MPN[i]*v[i])
      
      # Biological model for microbial abundance
      MPN[i] ~ dpois(lambda[i])
      log(lambda[i]) <- b0 + b1*watervib[water[i]] + U[site[i]] + V[samp[i]]
      
      # likelihood function for loo
      like[i] <- log(pbin(c[i],p[i],3))
      
    }
    for (s in 1:4) {
      U[s] ~ dnorm(0,tau_U)
    }
    for (t in 1:222) {
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
    for (i in 1:1362) {
      # Observation model across serial dilutions
      c[i] ~ dbin(p[i],3)
      p[i] <- 1-exp(-MPN[i]*v[i])
      
      # Biological model for microbial abundance
      MPN[i] ~ dpois(lambda[i])
      log(lambda[i]) <- b0 + b1*watervib[water[i]] + b2*temp[i] + U[site[i]] + V[samp[i]]
      
      # likelihood function for loo
      like[i] <- log(pbin(c[i],p[i],3))
      
    }
    for (s in 1:4) {
      U[s] ~ dnorm(0,tau_U)
    }
    for (t in 1:222) {
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
    for (i in 1:1362) {
      # Observation model across serial dilutions
      c[i] ~ dbin(p[i],3)
      p[i] <- 1-exp(-MPN[i]*v[i])
      
      # Biological model for microbial abundance
      MPN[i] ~ dpois(lambda[i])
      log(lambda[i]) <- b0 + b1*watervib[water[i]] + b2*temp[i] + b3*watervib[water[i]]*watervib[water[i]] + U[site[i]] + V[samp[i]]
      
      # likelihood function for loo
      like[i] <- log(pbin(c[i],p[i],3))
      
    }
    for (s in 1:4) {
      U[s] ~ dnorm(0,tau_U)
    }
    for (t in 1:222) {
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
inits.m0 <- list(list("U"=numeric(4),"V"=numeric(222),"tau_U"=0.1,"tau_V"=0.1,"b0"=1),
                 list("U"=numeric(4),"V"=numeric(222),"tau_U"=0.01,"tau_V"=0.1,"b0"=1),
                 list("U"=numeric(4),"V"=numeric(222),"tau_U"=1,"tau_V"=0.1,"b0"=1))
inits.m1 <- list(list("U"=numeric(4),"V"=numeric(222),"tau_U"=0.1,"tau_V"=0.1,"b0"=1,"b1"=0),
                 list("U"=numeric(4),"V"=numeric(222),"tau_U"=0.01,"tau_V"=0.1,"b0"=1,"b1"=0),
                 list("U"=numeric(4),"V"=numeric(222),"tau_U"=1,"tau_V"=0.1,"b0"=1,"b1"=0))
inits.m2 <- list(list("U"=numeric(4),"V"=numeric(222),"tau_U"=0.1,"tau_V"=0.1,"b0"=1,"b1"=0,"b2"=0),
                 list("U"=numeric(4),"V"=numeric(222),"tau_U"=0.01,"tau_V"=0.1,"b0"=1,"b1"=0,"b2"=0),
                 list("U"=numeric(4),"V"=numeric(222),"tau_U"=1,"tau_V"=0.1,"b0"=1,"b1"=0,"b2"=0))
inits.m3 <- list(list("U"=numeric(4),"V"=numeric(222),"tau_U"=0.1,"tau_V"=0.1,"b0"=1,"b1"=0,"b2"=0,"b3"=0),
                 list("U"=numeric(4),"V"=numeric(222),"tau_U"=0.01,"tau_V"=0.1,"b0"=1,"b1"=0,"b2"=0,"b3"=0),
                 list("U"=numeric(4),"V"=numeric(222),"tau_U"=1,"tau_V"=0.1,"b0"=1,"b1"=0,"b2"=0,"b3"=0))

params.base <- c("b0")
params.m0 <- c("b0")
params.m1 <- c("b0","b1")
params.m2 <- c("b0","b1","b2")
params.m3 <- c("b0","b1","b2","b3")

m <- jags(data = vibrio,
          inits = inits.m3,
          parameters.to.save = params.m3,
          model.file = "m3.jag",
          n.chains = 3,
          n.iter = 10000,
          n.burnin = 1000,
          n.thin = 3)

m.parmlist <- m$BUGSoutput$sims.list

# LOO
#m.loglike <- m.parmlist$log.like
#m.loo <- loo(m.loglike, r_eff = NA)
#m.loo
#m.mcmc <- as.mcmc(m)
#m.gel <- gelman.diag(m.mcmc, confidence = 0.95, transform=FALSE, autoburnin=TRUE,
#                     multivariate=TRUE)
#m.gel

# Parameter estimates and P
m.b0 <- m.parmlist$b0 
m.b0.P <- 1 - length(m.b0[m.b0>0])/length(m.b0)

m.b1 <- m.parmlist$b1 
m.b1.P <- 1 - length(m.b1[m.b1>0])/length(m.b1)

m.b2 <- m.parmlist$b2 
m.b2.P <- 1 - length(m.b2[m.b2>0])/length(m.b2)

m.b3 <- m.parmlist$b3
m.b3.P <- 1 - length(m.b3[m.b3>0])/length(m.b3)
m