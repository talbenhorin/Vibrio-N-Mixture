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

dat <- read.csv("seagrantvibrio.csv", fill = FALSE, header = TRUE) 
m1.dat<-list(c=dat$path,v=dat$mass,samp=dat$samp,gear=dat$gear,tide=dat$t2,mod=dat$mod,hi=dat$hi,time=dat$time)

# Mixed-effects model
# Only intercepts are random, but slopes are identical for all groups 

# N-Mixture model for serial dilution data 
cat(
  "model{
    for (i in 1:5320) {
      # Observation model across serial dilutions
      c[i] ~ dbin(p[i],3)
      p[i] <- 1-exp(-lambda[i]*v[i])
      
      # Biological model for microbial abundance
      lambda[i] ~ dlnorm(mu[i],tau_all)
      mu[i] <- b0 + b1*gear[i] + b2*tide[i] + b3*gear[i]*tide[i] + b4*mod[i] + b5*hi[i] + U[samp[i]] + V[time[i]]
    }
    for (s in 1:996) {
      U[s] ~ dnorm(0,tau_U)
    }
    for (t in 1:20) {
      V[t] ~ dnorm(0,tau_V)
    }
    tau_U ~ dgamma(0.1,0.1)
    tau_V ~ dgamma(0.1,0.1)
    tau_all ~ dgamma(0.1,0.1)
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
m1.inits <- list(list("U"=numeric(996),"V"=numeric(20),"tau_all"=0.1,"tau_U"=1,"tau_V"=0.1,"b0"=0,"b1"=0,"b2"=0,"b3"=0,"b4"=0,"b5"=0),
                 list("U"=numeric(996),"V"=numeric(20),"tau_all"=0.01,"tau_U"=1,"tau_V"=0.1,"b0"=0,"b1"=0,"b2"=0,"b3"=0,"b4"=0,"b5"=0),
                 list("U"=numeric(996),"V"=numeric(20),"tau_all"=0.01,"tau_U"=1,"tau_V"=0.1,"b0"=0,"b1"=0,"b2"=0,"b3"=0,"b4"=0,"b5"=0))

parameters <- c("U","V","b0","b1","b2","b3","b4","b5")

m1 <- jags(data = m1.dat,
           inits = m1.inits,
           parameters.to.save = parameters,
           model.file = "m1.jag",
           n.chains = 3,
           n.iter = 5000,
           n.burnin = 2000,
           n.thin = 3)  

# Posterior Median and Highest Posterior Density Intervals
out<-MCMCpstr(m1,
              params = parameters,
              func = median,
              type = 'summary')
out95<-hdi(list(m1$BUGSoutput$sims.list$b0,m1$BUGSoutput$sims.list$b1,m1$BUGSoutput$sims.list$b2,m1$BUGSoutput$sims.list$b3,m1$BUGSoutput$sims.list$b4,m1$BUGSoutput$sims.list$b5))
# output file
med <- rbind(out[3],out[4],out[5],out[6],out[7],out[8])
lower <- rbind(out95[[1]][1,],out95[[2]][1,],out95[[3]][1,],out95[[4]][1,],out95[[5]][1,],out95[[6]][1,])
upper <- rbind(out95[[1]][2,],out95[[2]][2,],out95[[3]][2,],out95[[4]][2,],out95[[5]][2,],out95[[6]][2,])
Iwant <-data.frame(as.numeric(med), lower, upper)
write.csv(Iwant,file = "output.csv", row.names = FALSE)
