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
library(ggplot2)

dat <- read.csv("seagrantvibrio.csv", fill = FALSE, header = TRUE) 

# N-Mixture model for serial dilution data 
cat(
  "model{
    for (i in 1:5320) {
      # Observation model across serial dilutions
      c[i] ~ dbin(p[i],3)
      p[i] <- 1-exp(-MPN[i]*v[i])
      
      # Biological model for microbial abundance
      MPN[i] ~ dpois(lambda[i])
      log(lambda[i]) <- b0 + b1*gear[i] + b2*tide[i] + b3*gear[i]*tide[i] + b4*mod[i] + b5*hi[i] + U[samp[i]] +V[time[i]]
    }
    for (s in 1:996) {
      U[s] ~ dnorm(0,tau_U)
    }
    for (t in 1:20) {
      V[t] ~ dnorm(0,tau_V)
    }
    tau_U ~ dgamma(0.1,0.1)
    tau_V ~ dgamma(0.1,0.1)
    b0 ~ dnorm(0,0.1)
    b1 ~ dnorm(0,0.1)
    b2 ~ dnorm(0,0.1)
    b3 ~ dnorm(0,0.1)
    b4 ~ dnorm(0,0.1)
    b5 ~ dnorm(0,0.1)
  }",
  file="m1.jag"
)

# Initial params BOTH YEARS
m1.inits <- list(list("U"=numeric(996),"V"=numeric(20),"tau_U"=0.1,"tau_V"=0.1,"b0"=0,"b1"=0,"b2"=0,"b3"=0,"b4"=0,"b5"=0),
                 list("U"=numeric(996),"V"=numeric(20),"tau_U"=0.01,"tau_V"=0.1,"b0"=0,"b1"=0,"b2"=0,"b3"=0,"b4"=0,"b5"=0),
                 list("U"=numeric(996),"V"=numeric(20),"tau_U"=1,"tau_V"=0.1,"b0"=0,"b1"=0,"b2"=0,"b3"=0,"b4"=0,"b5"=0))

parameters <- c("b0","b1","b2","b3","b4","b5")

m1.total <- list(c=dat$tlh,v=dat$mass,samp=dat$samp,gear=dat$gear,tide=dat$t2,mod=dat$mod,hi=dat$hi,time=dat$time) #data string, total vibrio
m1.path <-list(c=dat$path,v=dat$mass,samp=dat$samp,gear=dat$gear,tide=dat$t2,mod=dat$mod,hi=dat$hi,time=dat$time) #data string, pathogenic vibrio

mTotal <- jags(data = m1.total,
                inits = m1.inits,
                parameters.to.save = parameters,
                model.file = "m1.jag",
                n.chains = 3,
                n.iter = 4000,
                n.burnin = 1000,
                n.thin = 3)  

mPath <- jags(data = m1.path,
                inits = m1.inits,
                parameters.to.save = parameters,
                model.file = "m1.jag",
                n.chains = 3,
                n.iter = 4000,
                n.burnin = 1000,
                n.thin = 3)

out <- data.frame(
  Total=mTotal$BUGSoutput$sims.list$b3,
  Pathogenic=mPath$BUGSoutput$sims.list$b3)
b3 <- stack(out,select=c('Total','Pathogenic'))
names(b3) <- c("values","vibrio")

ggplot(b3, aes(x=values, fill=vibrio))+ 
  geom_density(alpha=0.4)+
  theme_light()+
  labs(title=expression(beta[3]),x = "", y = "")
  



  

