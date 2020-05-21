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

dat <- read.csv("va2014.csv", fill = FALSE, header = TRUE) 

# Base and Full N-Mixture model for serial dilution data 
cat(
  "model{
    for (i in 1:288) {
      # Observation model across serial dilutions
      c[i] ~ dbin(p[i],3)
      p[i] <- 1-exp(-MPN[i]*v[i])
      
      # Biological model for microbial abundance
      MPN[i] ~ dpois(lambda[i])
      log(lambda[i]) <- b0 + U[site[i]] + V[samp[i]]
      log.like[i] <- log(pbin(c[i],p[i],3))*log(ppois(MPN[i],lambda[i]))
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
      log(lambda[i]) <- b0 + b1*temp[i] + b2*chlo[i] + b3*temp[i]*chlo[i] + U[site[i]] + V[samp[i]]
      log.like[i] <- log(pbin(c[i],p[i],3))*log(ppois(MPN[i],lambda[i]))
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
  file="full.jag"
)

# Initial params BOTH YEARS
base.inits <- list(list("U"=numeric(3),"V"=numeric(48),"tau_U"=0.1,"tau_V"=0.1,"b0"=0),
                   list("U"=numeric(3),"V"=numeric(48),"tau_U"=0.01,"tau_V"=0.1,"b0"=0),
                   list("U"=numeric(3),"V"=numeric(48),"tau_U"=1,"tau_V"=0.1,"b0"=0))


full.inits <- list(list("U"=numeric(3),"V"=numeric(48),"tau_U"=0.1,"tau_V"=0.1,"b0"=0,"b1"=0,"b2"=0,"b3"=0),
                 list("U"=numeric(3),"V"=numeric(48),"tau_U"=0.01,"tau_V"=0.1,"b0"=0,"b1"=0,"b2"=0,"b3"=0),
                 list("U"=numeric(3),"V"=numeric(48),"tau_U"=1,"tau_V"=0.1,"b0"=0,"b1"=0,"b2"=0,"b3"=0))

parameters <- c("log.like")
pbase <- c("b0")
pfull <- c("b0","b1","b2","b3")

Vv.total <- list(c=dat$Vv.vvha,v=dat$Sample.Volume,samp=dat$FID,site=dat$Site.Num,temp=dat$temp.t,pheo=dat$pheo.t,turb=dat$turb.t,chlo=dat$chlo.t) #data string, total vibrio
Vv.path <-list(c=dat$Vv.PilF,v=dat$Sample.Volume,samp=dat$FID,site=dat$Site.Num,temp=dat$temp.t,pheo=dat$pheo.t,turb=dat$turb.t,chlo=dat$chlo.t) #data string, pathogenic vibrio
Vp.total <- list(c=dat$Vp.tlh,v=dat$Sample.Volume,samp=dat$FID,site=dat$Site.Num,temp=dat$temp.t,pheo=dat$pheo.t,turb=dat$turb.t,chlo=dat$chlo.t) #data string, total vibrio
Vp.path <- list(c=dat$Vp.path,v=dat$Sample.Volume,samp=dat$FID,site=dat$Site.Num,temp=dat$temp.t,pheo=dat$pheo.t,turb=dat$turb.t,chlo=dat$chlo.t) #data string, total vibrio

m.base <- jags(data = Vp.path,
               inits = base.inits,
               parameters.to.save = parameters,
               model.file = "base.jag",
               n.chains = 3,
               n.iter = 10000,
               n.burnin = 1000,
               n.thin = 3)
m.full <- jags(data = Vp.total,
               inits = full.inits,
               parameters.to.save = pfull,
               model.file = "full.jag",
               n.chains = 3,
               n.iter = 10000,
               n.burnin = 1000,
               n.thin = 3)

# WAIC
base.waic <- waic(m.base$BUGSoutput$sims.list$log.like)
full.waic <- waic(m.full$BUGSoutput$sims.list$log.like)
full.waic
compare(base.waic, full.waic)
#write.csv(ln.waic$estimates,file = "WAIC.csv", row.names = FALSE)

