rm(list=ls(all=TRUE))

library(R2jags)
library(loo)
library(HDInterval)
library(MCMCvis)

dat <- read.csv("yr2seagrantvibrio.csv", fill = FALSE, header = TRUE) 
m1.dat<-list(c=dat$tlh,v=dat$mass,samp=dat$samp,gear=dat$gear,tide=dat$tmax)

# N-Mixture model for serial dilution data 
cat(
  "model{
    for (i in 1:920) {
      # Observation model across serial dilutions
      c[i] ~ dbin(p[i],3)
      p[i] <- 1-exp(-MPN[i]*v[i])
      
      # Biological model for microbial abundance
      MPN[i] ~ dpois(lambda[i])
      log(lambda[i]) <- b0 + b1*gear[i] + b2*tide[i] + b3*gear[i]*tide[i] + U[samp[i]] 
    }
    for (s in 1:168) {
      U[s] ~ dnorm(0,tau_U)
    }
    tau_U ~ dgamma(0.1,0.1)
    b0 ~ dnorm(0,0.1)
    b1 ~ dnorm(0,0.1)
    b2 ~ dnorm(0,0.1)
    b3 ~ dnorm(0,0.1)
  }",
  file="m1.jag"
)

# Initial params BOTH YEARS
m1.inits <- list(list("U"=numeric(168),"tau_U"=0.1,"b0"=0,"b1"=0,"b2"=0,"b3"=0),
                 list("U"=numeric(168),"tau_U"=0.01,"b0"=0,"b1"=0,"b2"=0,"b3"=0),
                 list("U"=numeric(168),"tau_U"=1,"b0"=0,"b1"=0,"b2"=0,"b3"=0))

parameters <- c("U","b0","b1","b2","b3")

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
out95<-hdi(list(m1$BUGSoutput$sims.list$b0,
                m1$BUGSoutput$sims.list$b1,
                m1$BUGSoutput$sims.list$b2,
                m1$BUGSoutput$sims.list$b3))

# output file
med <- rbind(out[2],out[3],out[4],out[5])
lower <- rbind(out95[[1]][1,],
               out95[[2]][1,],
               out95[[3]][1,],
               out95[[4]][1,])
upper <- rbind(out95[[1]][2,],
               out95[[2]][2,],
               out95[[3]][2,],
               out95[[4]][2,])
Iwant <-data.frame(as.numeric(med), lower, upper)   
write.csv(Iwant,file = "output.csv", row.names = FALSE)