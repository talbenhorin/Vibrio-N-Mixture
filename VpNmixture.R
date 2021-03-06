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
    LOR0 <- exp(b0)
    b1 ~ dnorm(0,0.1)
    LOR1 <- exp(b1)
    b2 ~ dnorm(0,0.1)
    LOR2 <- exp(b2)
    b3 ~ dnorm(0,0.1)
    LOR3 <- exp(b3)
    b4 ~ dnorm(0,0.1)
    LOR4 <- exp(b4)
    b5 ~ dnorm(0,0.1)
    LOR5 <- exp(b5)
  }",
  file="m1.jag"
)

# Initial params BOTH YEARS
m1.inits <- list(list("U"=numeric(996),"V"=numeric(20),"tau_U"=0.1,"tau_V"=0.1,"b0"=0,"b1"=0,"b2"=0,"b3"=0,"b4"=0,"b5"=0),
                 list("U"=numeric(996),"V"=numeric(20),"tau_U"=0.01,"tau_V"=0.1,"b0"=0,"b1"=0,"b2"=0,"b3"=0,"b4"=0,"b5"=0),
                 list("U"=numeric(996),"V"=numeric(20),"tau_U"=1,"tau_V"=0.1,"b0"=0,"b1"=0,"b2"=0,"b3"=0,"b4"=0,"b5"=0))

parameters <- c("LOR0","LOR1","LOR2","LOR3","LOR4","LOR5")

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

out0 <- data.frame(
  Total=mTotal$BUGSoutput$sims.list$LOR0,
  Pathogenic=mPath$BUGSoutput$sims.list$LOR0)
b0 <- stack(out0,select=c('Total','Pathogenic'))
names(b0) <- c("values","vibrio")

out1 <- data.frame(
  Total=mTotal$BUGSoutput$sims.list$LOR1,
  Pathogenic=mPath$BUGSoutput$sims.list$LOR1)
b1 <- stack(out1,select=c('Total','Pathogenic'))
names(b1) <- c("values","vibrio")

out2 <- data.frame(
  Total=mTotal$BUGSoutput$sims.list$LOR2,
  Pathogenic=mPath$BUGSoutput$sims.list$LOR2)
b2 <- stack(out2,select=c('Total','Pathogenic'))
names(b2) <- c("values","vibrio")

out3 <- data.frame(
  Total=mTotal$BUGSoutput$sims.list$LOR3,
  Pathogenic=mPath$BUGSoutput$sims.list$LOR3)
b3 <- stack(out3,select=c('Total','Pathogenic'))
names(b3) <- c("values","vibrio")

out4 <- data.frame(
  Total=mTotal$BUGSoutput$sims.list$LOR4,
  Pathogenic=mPath$BUGSoutput$sims.list$LOR4)
b4 <- stack(out4,select=c('Total','Pathogenic'))
names(b4) <- c("values","vibrio")

out5 <- data.frame(
  Total=mTotal$BUGSoutput$sims.list$LOR5,
  Pathogenic=mPath$BUGSoutput$sims.list$LOR5)
b5 <- stack(out5,select=c('Total','Pathogenic'))
names(b5) <- c("values","vibrio")

p0 <- ggplot(b0, aes(x=values, fill=vibrio))+ 
  geom_density(alpha=0.4)+
  theme_classic()+
  scale_fill_grey()+
  scale_x_log10(breaks = c(1,10,100,1000),
                labels = c("1" = "1","10" = "10","100"="100","1000"="1000"))+
  labs(x = expression(CFU~g^-1), y = "")+
  theme(text=element_text(size=14,  family="Helvetica"))+
  theme(legend.justification=c(1,0), legend.position=c(0.19,0.65))+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank())

p1 <- ggplot(b1, aes(x=values, fill=vibrio))+ 
  geom_density(alpha=0.4)+
  theme_void()+
  scale_fill_grey()+
  theme(
    legend.position="none",strip.text.y=element_text(angle=0, hjust=0.5),
    panel.border=element_rect(fill=NA, color = "white", size=0.67),
    plot.margin=margin(t=2, r=4, b=2, l=2, unit="pt"),
    panel.grid=element_blank(), panel.ontop=FALSE, axis.title.y=element_blank())+
  xlim(0,3)+
  labs(x = "", y = "")+
  theme(plot.title = element_text(hjust = 0))

p2 <- ggplot(b2, aes(x=values, fill=vibrio))+ 
  geom_density(alpha=0.4)+
  theme_void()+
  scale_fill_grey()+
  theme(
    legend.position="none",strip.text.y=element_text(angle=0, hjust=0.5),
    panel.border=element_rect(fill=NA, color = "white", size=0.67),
    plot.margin=margin(t=2, r=4, b=2, l=2, unit="pt"),
    panel.grid=element_blank(), panel.ontop=FALSE, axis.title.y=element_blank())+
  xlim(0,3)+
  labs(x = "", y = "")+
  theme(plot.title = element_text(hjust = 0))

p3 <- ggplot(b3, aes(x=values, fill=vibrio))+ 
  geom_density(alpha=0.4)+
  theme_void()+
  scale_fill_grey()+
  theme(
    legend.position="none",strip.text.y=element_text(angle=0, hjust=0.5),
    panel.border=element_rect(fill=NA, color = "white", size=0.67),
    plot.margin=margin(t=2, r=4, b=2, l=2, unit="pt"),
    panel.grid=element_blank(), panel.ontop=FALSE, axis.title.y=element_blank())+
  xlim(0,3)+
  labs(x = "", y = "")+
  theme(plot.title = element_text(hjust = 0))

p4 <- ggplot(b4, aes(x=values, fill=vibrio))+ 
  geom_density(alpha=0.4)+
  theme_void()+
  scale_fill_grey()+
  theme(
    legend.position="none",strip.text.y=element_text(angle=0, hjust=0.5),
    panel.border=element_rect(fill=NA, color = "white", size=0.67),
    plot.margin=margin(t=2, r=4, b=2, l=2, unit="pt"),
    panel.grid=element_blank(), panel.ontop=FALSE, axis.title.y=element_blank())+
  xlim(0,3)+
  labs(x = "", y = "")+
  theme(plot.title = element_text(hjust = 0))

  p5 <- ggplot(b5, aes(x=values, fill=vibrio))+ 
  geom_density(alpha=0.4)+
  theme_classic()+
  scale_fill_grey()+
  theme(
    text=element_text(size=14,  family="Helvetica"),
    legend.position="none",strip.text.y=element_text(angle=0, hjust=0.5),
    panel.border=element_rect(fill=NA, color = "white", size=0.67),
    plot.margin=margin(t=2, r=4, b=2, l=2, unit="pt"),
    panel.grid=element_blank(), panel.ontop=FALSE, axis.title.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank())+
  theme(plot.title = element_text(hjust = 0))+
  xlim(0,3)+
  labs(x = expression("Relative Risk"~(e^beta)), y = "")

(p0 / p1 / p2 / p3 / p4 / p5)



