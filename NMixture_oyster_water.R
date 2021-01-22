rm(list=ls(all=TRUE))

# Download JAGS at https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Windows/
# Load these packages the first time your run R
install.packages("RTools")
install.packages("rjags")
install.packages("coda")
install.packages("R2jags")
install.packages("hdi")
install.packages("MCMCvis")
install.packages("loo")
install.packages("HDInterval")

library(R2jags)
library(loo)
library(HDInterval)
library(MCMCvis)
library(ggplot2)
library(patchwork)
library(scales)

dat <- read.csv("va2014water.csv", fill = FALSE, header = TRUE) 

# N-Mixture model for serial dilution data 
cat(
  "model{
<<<<<<< HEAD
    for (i in 1:288) {
=======
    for (i in 1:1362) {
>>>>>>> bb652fa6ef87a040c32aeaf636efacea854e4049
      # Observation model across serial dilutions
      c[i] ~ dbin(p[i],3)
      p[i] <- 1-exp(-MPN[i]*v[i])
      
      # Biological model for microbial abundance
      MPN[i] ~ dpois(lambda[i])
<<<<<<< HEAD
      log(lambda[i]) <- b0 + b1*W[water[i]] + U[site[i]] + V[samp[i]]
      post.prob[i] <- pbin(c[i],p[i],3) + ppois(MPN[i],lambda[i])
=======
      log(lambda[i]) <- b0 + b1*chlo[i] + U[site[i]] + V[samp[i]]
      log.like[i] <- log(pbin(c[i],p[i],3))*log(ppois(MPN[i],lambda[i]))
>>>>>>> edd59bb33139fc9a512b0c84dda658927cb111ea
    }
    for (s in 1:3) {
      U[s] ~ dnorm(0,tau_U)
    }
<<<<<<< HEAD
    for (t in 1:48) {
      V[t] ~ dnorm(0,tau_V)
    }
=======
    for (t in 1:222) {
      V[t] ~ dnorm(0,tau_V)
    }
    for (w in 1:56) {
      W[w] <- mu_w[w]
    }
>>>>>>> bb652fa6ef87a040c32aeaf636efacea854e4049
    tau_U ~ dgamma(0.1,0.1)
    tau_V ~ dgamma(0.1,0.1)
    b0 ~ dnorm(0,0.1)
    b1 ~ dnorm(0,0.1)
    r ~ dunif(0,50)
  }",
  file="water.jag"
)

# Initial params BOTH YEARS
<<<<<<< HEAD
inits <- list(list("U"=numeric(3),"V"=numeric(48),"tau_U"=0.1,"tau_V"=0.1,"b0"=0,"b1"=0),
                   list("U"=numeric(3),"V"=numeric(48),"tau_U"=0.01,"tau_V"=0.1,"b0"=0,"b1"=0),
                   list("U"=numeric(3),"V"=numeric(48),"tau_U"=1,"tau_V"=0.1,"b0"=0,"b1"=0))

pfull <- c("b0","b1","log.like")

in.data <- list(c=dat$Vv.vvha,v=dat$Sample.Volume,samp=dat$FID,site=dat$Site.Num,chlo=dat$chlo.t) #data string
=======
inits <- list(list("U"=numeric(4),"V"=numeric(222),"tau_U"=0.1,"tau_V"=0.1,"b0"=0,"b1"=0),
                   list("U"=numeric(4),"V"=numeric(222),"tau_U"=0.01,"tau_V"=0.1,"b0"=0,"b1"=0),
                   list("U"=numeric(4),"V"=numeric(222),"tau_U"=1,"tau_V"=0.1,"b0"=0,"b1"=0))

pfull <- c("post.prob")

in.data <- list(c=dat$pilf,v=dat$mass,samp=dat$fid,site=dat$site,water=dat$water,mu_w=dat$water.vvha) #data string
>>>>>>> bb652fa6ef87a040c32aeaf636efacea854e4049

m.base <- jags(data = in.data,
               inits = inits,
               parameters.to.save = pfull,
               model.file = "water.jag",
               n.chains = 3,
               n.iter = 4000,
               n.burnin = 1000,
               n.thin = 3)


<<<<<<< HEAD


=======
#out <-MCMCpstr(m.base,
#              params = pfull,
#              func = median,
#              type = 'summary')
b0.95 <- hdi(m.base$BUGSoutput$sims.list$b0)
b1.95 <- hdi(m.base$BUGSoutput$sims.list$b1)
<<<<<<< HEAD

=======
>>>>>>> bb652fa6ef87a040c32aeaf636efacea854e4049
>>>>>>> edd59bb33139fc9a512b0c84dda658927cb111ea

