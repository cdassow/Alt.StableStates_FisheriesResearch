## CD 12.22.2020
## Supplementary information figures
rm(list=ls())
library(deSolve)
library(ggplot2)
library(ggpubr)

#model function
# qE - harvest, species specific, this is what is controlled by 'regulations'
# s - juvenile overwinter survival, species specific
# m - adult natural mortality rate, species specific
# cJA - effect of adults of a given species on juveniles of a given species (cover cannibalism or interspecific predation, both happen in foraging arena)
# cJJ - effect of juveniles of one species on juveniles of the other (can be predation or competition)
# h - rate at which juveniles leave foraging arena for refuge, species specific
# v - rate at which juveniles enter foraging arena from refuge, species specific
# stock1 - annual stocked num spp1
# stock2 - annual stocked num spp2
simBiggsQ2<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
    dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
    dJ1dt=-cJ1J2*J2*J1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))-((cJ1A1*v1*A1*J1)/(h1Fun(t)+v1+cJ1A1*A1))-s1*J1+(100*A1/(50+A1))+st1Fun(t)
    dJ2dt=-cJ2J1*J1*J2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))-((cJ2A2*v2*A2*J2)/(h2Fun(t)+v2+cJ2A2*A2))-s2*J2+(100*A2/(50+A2))+st2Fun(t)
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}



 ## SI figure 2
#flip scenarios

tstep=1:100

h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
minDiff=100

#### FLIP W/O sp2 harv ####

qEs=seq(0,0.3,length.out = 30)
sto=seq(0,200, length.out = 30)
df2=expand.grid(X=qEs,Y=sto)
df2$A1=numeric(nrow(df2))
df2$A2=numeric(nrow(df2))
df2$J1=numeric(nrow(df2))
df2$J2=numeric(nrow(df2))


for(i in 1:nrow(df2)){
  tstep=1:100
  qE1Fun=approxfun(x=tstep,y=rep(df2$X[i], length(tstep)))
  qE2Fun=approxfun(x=tstep,y=rep(0, length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(df2$Y[i],length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  p=c(s1=0.1,m1=0.1,cJ1A1=0.002,cJ1A2=0.05,cJ1J2=0.003,v1=1,
      s2=0.1,m2=0.1,cJ2A2=0.002,cJ2A1=0.03,cJ2J1=0.003,v2=1)
  y0=c(500,5000,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  df2$A1[i]=sim[nrow(sim)-1,2]
  df2$A2[i]=sim[nrow(sim)-1,3]
  df2$J1[i]=sim[nrow(sim)-1,4]
  df2$J2[i]=sim[nrow(sim)-1,5]
}

df2$diff=df2$A1-df2$A2
df2$ratio=df2$A2/df2$A1
df2$outcome=ifelse(df2$ratio < 0.6,"darkgreen", "darkred")


#### FLIP W/ sp2 harv ####

qEs=seq(0,0.3,length.out = 30)
sto=seq(0,200, length.out = 30)
dfwo2=expand.grid(X=qEs,Y=sto)
dfwo2$A1=numeric(nrow(dfwo2))
dfwo2$A2=numeric(nrow(dfwo2))
dfwo2$J1=numeric(nrow(dfwo2))
dfwo2$J2=numeric(nrow(dfwo2))

for(i in 1:nrow(dfwo2)){
  tstep=1:100
  qE1Fun=approxfun(x=tstep,y=rep(dfwo2$X[i], length(tstep)))
  qE2Fun=approxfun(x=tstep,y=rep(0.1, length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(dfwo2$Y[i],length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  p=c(s1=0.1,m1=0.1,cJ1A1=0.002,cJ1A2=0.05,cJ1J2=0.003,v1=1,
      s2=0.1,m2=0.1,cJ2A2=0.002,cJ2A1=0.03,cJ2J1=0.003,v2=1)
  y0=c(500,5000,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  dfwo2$A1[i]=sim[nrow(sim)-1,2]
  dfwo2$A2[i]=sim[nrow(sim)-1,3]
  dfwo2$J1[i]=sim[nrow(sim)-1,4]
  dfwo2$J2[i]=sim[nrow(sim)-1,5]
}

dfwo2$diff=dfwo2$A1-dfwo2$A2
dfwo2$ratio=dfwo2$A2/dfwo2$A1
dfwo2$outcome=ifelse(dfwo2$ratio < 0.6,"darkgreen", "darkred")
df2$mod=rep("Flip to Preferred, ignore 2",nrow(df2))
dfwo2$mod=rep("Flip to Preferred, harv 2",nrow(dfwo2))
allSen=rbind(df2,dfwo2)

## use this plot for SI figure 2
ggplot(data=allSen, aes(x=allSen$X,y=allSen$Y,linetype=allSen$mod))+theme_classic()+
  geom_contour(aes(z=allSen$ratio),breaks = c(0.6), color='black', size=1.7)+
  labs(x="Preferred Species Harvest Rate", y="Preferred Species Stocking",linetype="Scenario")+
  theme(legend.position = 'none')+
  geom_segment(aes(x=0.15, y=65, xend=0.2, yend=40), arrow = arrow(type='closed'),size=1)


## SI FIGURE 3
#plot to look at the cost/benefits of stocking or predator reduction for managing a focal species. Flip scenario instead of maintain scenario 

minDiff=100
qEs=rep(seq(0,0.3,length.out = 30),3)
sto=rep(seq(0,200, length.out = 30),3)
dfT=expand.grid(X=qEs,Y=sto)
dfT$A1=numeric(nrow(dfT))
dfT$A2=numeric(nrow(dfT))
dfT$J1=numeric(nrow(dfT))
dfT$J2=numeric(nrow(dfT))
dfT$sp1H=c(rep(0.05,30),rep(0.1,30),rep(0.15,30))

for(i in 1:nrow(dfT)){
  tstep=1:100
  qE1Fun=approxfun(x=tstep,y=rep(dfT$sp1H[i], length(tstep)))
  qE2Fun=approxfun(x=tstep,y=rep(dfT$X[i], length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(dfT$Y[i],length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  p=c(s1=0.1,m1=0.1,cJ1A1=0.002,cJ1A2=0.05,cJ1J2=0.003,v1=1,
      s2=0.1,m2=0.1,cJ2A2=0.002,cJ2A1=0.03,cJ2J1=0.003,v2=1)
  y0=c(500,5000,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  dfT$A1[i]=sim[nrow(sim)-1,2]
  dfT$A2[i]=sim[nrow(sim)-1,3]
  dfT$J1[i]=sim[nrow(sim)-1,4]
  dfT$J2[i]=sim[nrow(sim)-1,5]
}

dfT$diff=dfT$A1-dfT$A2
dfT$ratio=dfT$A2/dfT$A1
dfT$qNorm=(dfT$X-min(dfT$X))/(max(dfT$X)-min(dfT$X))
dfT$sNorm=(dfT$Y-min(dfT$Y))/(max(dfT$Y)-min(dfT$Y))
dfT$sp1Norm=((dfT$sp1H-min(dfT$X))/(max(dfT$X)-min(dfT$X))) #putting sp1 harv on same scale as sp2

# SI figure 3
vzT=ggplot(data=dfT, aes(x=dfT$X,y=dfT$Y,linetype=as.factor(dfT$sp1H)))+theme_classic()+
  geom_contour(aes(z=dfT$ratio),breaks = c(0.6), color='black', size=1)+
  labs(x="Competitor Species Harvest Rate", y="Preferred Species Stocking",linetype="Preferred Species Harvest")+
  theme(legend.position = 'bottom')
vzT
