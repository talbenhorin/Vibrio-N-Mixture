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

dat <- read.csv("vaoyster.csv", fill = FALSE, header = TRUE) 

# N-Mixture model for serial dilution data 
cat(
  "model{
    for (i in 1:1796) {
      # Observation model across serial dilutions
      c[i] ~ dbin(p[i],3)
      p[i] <- 1-exp(-MPN[i]*v[i])
      
      # Biological model for microbial abundance
      MPN[i] ~ dpois(lambda[i])
      log(lambda[i]) <- b0 + b1*W[water[i]] + U[site[i]] + V[samp[i]]
      log.like[i] <- log(pbin(c[i],p[i],3))*log(ppois(MPN[i],lambda[i]))
    }
    for (s in 1:4) {
      U[s] ~ dnorm(0,tau_U)
    }
    for (t in 1:284) {
      V[t] ~ dnorm(0,tau_V)
    }
    for (w in 1:71) {
      W[w] <- mu_w[w]
    }
    tau_U ~ dgamma(0.1,0.1)
    tau_V ~ dgamma(0.1,0.1)
    b0 ~ dnorm(0,0.1)
    b1 ~ dnorm(0,0.1)
  }",
  file="water.jag"
)

# Initial params BOTH YEARS
inits <- list(list("U"=numeric(4),"V"=numeric(284),"tau_U"=0.1,"tau_V"=0.1,"b0"=0,"b1"=-0),
                   list("U"=numeric(4),"V"=numeric(284),"tau_U"=0.01,"tau_V"=0.1,"b0"=0,"b1"=0),
                   list("U"=numeric(4),"V"=numeric(284),"tau_U"=1,"tau_V"=0.1,"b0"=0,"b1"=0))

pfull <- c("b0","b1","log.like")

in.data <- list(c=dat$vvha,v=dat$mass,samp=dat$fid,site=dat$site,water=dat$water,mu_w=dat$water.vvha) #data string

m.base <- jags(data = in.data,
               inits = inits,
               parameters.to.save = pfull,
               model.file = "water.jag",
               n.chains = 3,
               n.iter = 10000,
               n.burnin = 1000,
               n.thin = 3)

m.waic <- waic(m.base$BUGSoutput$sims.list$log.like)

out <-MCMCpstr(m.base,
              params = pfull,
              func = median,
              type = 'summary')
b0.95 <- hdi(m.base$BUGSoutput$sims.list$b0)
b1.95 <- hdi(m.base$BUGSoutput$sims.list$b1)