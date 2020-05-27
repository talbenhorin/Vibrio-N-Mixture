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

dall <- read.csv("vaoyster.csv", fill = FALSE, header = TRUE) 
dpilf <- read.csv("vaoysterpilf.csv", fill = FALSE, header = TRUE) 

cat(
  "model{
    for (i in 1:1796) {
      # Observation model across serial dilutions
      c[i] ~ dbin(p[i],3)
      p[i] <- 1-exp(-MPN[i]*v[i])
      
      # Biological model for microbial abundance
      MPN[i] ~ dpois(lambda[i])
      log(lambda[i]) <- b0 + b1*W[water[i]] + U[site[i]] + V[samp[i]]
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

cat(
  "model{
    for (i in 1:1404) {
      # Observation model across serial dilutions
      c[i] ~ dbin(p[i],3)
      p[i] <- 1-exp(-MPN[i]*v[i])
      
      # Biological model for microbial abundance
      MPN[i] ~ dpois(lambda[i])
      log(lambda[i]) <- b0 + b1*W[water[i]] + U[site[i]] + V[samp[i]]
    }
    for (s in 1:4) {
      U[s] ~ dnorm(0,tau_U)
    }
    for (t in 1:228) {
      V[t] ~ dnorm(0,tau_V)
    }
    for (w in 1:57) {
      W[w] <- mu_w[w]
    }
    tau_U ~ dgamma(0.1,0.1)
    tau_V ~ dgamma(0.1,0.1)
    b0 ~ dnorm(0,0.1)
    b1 ~ dnorm(0,0.1)
  }",
  file="waterpilf.jag"
)

# Initial params BOTH YEARS
inits <- list(list("U"=numeric(4),"V"=numeric(284),"tau_U"=0.1,"tau_V"=0.1,"b0"=0,"b1"=0),
                  list("U"=numeric(4),"V"=numeric(284),"tau_U"=0.01,"tau_V"=0.1,"b0"=0,"b1"=0),
                  list("U"=numeric(4),"V"=numeric(284),"tau_U"=1,"tau_V"=0.1,"b0"=0,"b1"=0))

initspilf <- list(list("U"=numeric(4),"V"=numeric(228),"tau_U"=0.1,"tau_V"=0.1,"b0"=0,"b1"=0),
              list("U"=numeric(4),"V"=numeric(228),"tau_U"=0.01,"tau_V"=0.1,"b0"=0,"b1"=0),
              list("U"=numeric(4),"V"=numeric(228),"tau_U"=1,"tau_V"=0.1,"b0"=0,"b1"=0))

params <- c("b0","b1")

in1 <- list(c=dall$vvha,v=dall$mass,samp=dall$fid,site=dall$site,water=dall$water,mu_w=dall$water.vvha) #data string
in2 <- list(c=dall$pilf,v=dall$mass,samp=dall$fid,site=dall$site,water=dall$water,mu_w=dall$water.vvha)
in3 <- list(c=dpilf$pilf,v=dpilf$mass,samp=dpilf$fid.new,site=dpilf$site,water=dpilf$water,mu_w=dpilf$water.plif)  
  
m1 <- jags(data = in1,
               inits = inits,
               parameters.to.save = params,
               model.file = "water.jag",
               n.chains = 3,
               n.iter = 5000,
               n.burnin = 2000,
               n.thin = 3)  
m2 <- jags(data = in2,
           inits = inits,
           parameters.to.save = params,
           model.file = "water.jag",
           n.chains = 3,
           n.iter = 5000,
           n.burnin = 2000,
           n.thin = 3) 
m3 <- jags(data = in3,
           inits = initspilf,
           parameters.to.save = params,
           model.file = "waterpilf.jag",
           n.chains = 3,
           n.iter = 4000,
           n.burnin = 2000,
           n.thin = 3) 

# Sort output for plot
out1 <- data.frame(o1=m1$BUGSoutput$sims.list$b1)
mod1 <- stack(out1,select=c('o1'))
names(mod1) <- c("pars","group")

out2 <- data.frame(o2=m2$BUGSoutput$sims.list$b1)
mod2 <- stack(out2,select=c('o2'))
names(mod2) <- c("pars","group")

out3 <- data.frame(o3=m3$BUGSoutput$sims.list$b1)
mod3 <- stack(out3,select=c('o3'))
names(mod3) <- c("pars","group")

# Plots
p1 <- ggplot(mod1, aes(pars, fill = group))+
  geom_density(alpha=0.25)+
  theme_classic()+
  theme(legend.position="none")+
  xlim(-2,3)+
  labs(title=expression("Water total Vv on oyster total Vv"),x = "", y = "")+
  geom_vline(xintercept = 0, linetype="dotted", 
             color = "black", size=1)
p2 <- ggplot(mod2, aes(pars, fill = group))+
  geom_density(alpha=0.25)+
  theme_classic()+
  theme(legend.position="none")+
  xlim(-2,3)+
  labs(title=expression("Water total Vv on oyster pathogenic Vv"),x = "", y = "")+
  geom_vline(xintercept = 0, linetype="dotted", 
             color = "black", size=1)
p3 <- ggplot(mod3, aes(pars, fill = group))+
  geom_density(alpha=0.25)+
  theme_classic()+
  theme(legend.position="none")+
  xlim(-2,3)+
  labs(title=expression("Water pathogenic Vv on oyster pathogenic Vv"),x = "Effect size", y = "")+
  geom_vline(xintercept = 0, linetype="dotted", 
             color = "black", size=1)
(p1 / p2 / p3)
  