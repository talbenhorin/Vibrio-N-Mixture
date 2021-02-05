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

dat <- read.csv("vaAllwater.csv", fill = FALSE, header = TRUE) 
Vv.total <- list(c=dat$Vp.path,v=dat$Sample.Volume,samp=dat$FID,site=dat$Site.Num,temp=dat$temp.t,pheo=dat$pheo.t,turb=dat$turb.t,chlo=dat$chlo.t) #data string, total vibrio
Vv.path <-list(c=dat$Vv.PilF,v=dat$Sample.Volume,samp=dat$FID,site=dat$Site.Num,temp=dat$temp.t,pheo=dat$pheo.t,turb=dat$turb.t,chlo=dat$chlo.t) #data string, pathogenic vibrio

# log.like[i] <- log(pbin(c[i],p[i],3))*log(ppois(MPN[i],lambda[i]))

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
      
      # likelihood function for loo
      like[i] <- log(pbin(c[i],p[i],3))

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
    
    #sum log likihoods for loo
    log.like <- sum(like)
    
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


params <- c("b0","b1","b2","b3")

m <- jags(data = Vv.total,
               inits = base.inits,
               parameters.to.save = params,
               model.file = "model1.jag",
               n.chains = 3,
               n.iter = 10000,
               n.burnin = 1000,
               n.thin = 3)

m.parmlist <- m$BUGSoutput$sims.list
m.loglike <- m.parmlist$log.like 
dim(m.loglike)

# GR: Gelman Rubin convergence statistic
# ESS: Effective samples size of 1,000,000 MCMC iterations

out<-MCMCpstr(m,
              params = params,
              func = median,
              type = 'summary')
out95<-hdi(list(m$BUGSoutput$sims.list$b0,
                m$BUGSoutput$sims.list$b1,
                m$BUGSoutput$sims.list$b2,
                m$BUGSoutput$sims.list$b3))

# output file
med <- rbind(out[1],out[2],out[3],out[4])
lower <- rbind(out95[[1]][1,],
               out95[[2]][1,],
               out95[[3]][1,],
               out95[[4]][1,])
upper <- rbind(out95[[1]][2,],
               out95[[2]][2,],
               out95[[3]][2,],
               out95[[4]][2,])
Iwant <-data.frame(as.numeric(med), lower, upper)   
Iwant
