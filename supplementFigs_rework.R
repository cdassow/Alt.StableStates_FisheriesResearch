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
# v. - rate at which juveniles leave foraging arena for refuge, species specific (equivalent to v' in foraging arena theory but easier to write in code because ' has a specific use). In the manuscript I will use the traditional v' notation.
# v - rate at which juveniles enter foraging arena from refuge, species specific
# kP - annual stocked num preferred species
# kC - annual stocked num competitor species
# wP - annual stochasticity added to natural recruitment for preferred species
# wC - annual stochasticity added to natural recruitment for competitor species
simBiggsQ2w<-function(t,y,params){
  AP<-y[1]
  AC<-y[2]
  JP<-y[3]
  JC<-y[4]
  with(as.list(params),{
    dAPdt=-qEPFun(t)*AP-mP*AP+sP*JP
    dACdt=-qECFun(t)*AC-mC*AC+sC*JC
    dJPdt=-cJPJC*JC*JP-((cJPAC*vP*AC*JP)/(v.PFun(t)+vP+cJPAC*AC))-((cJPAP*vP*AP*JP)/(v.PFun(t)+vP+cJPAP*AP))-sP*JP+((100*AP/(50+AP))*wPFun(t))+kPFun(t)
    dJCdt=-cJCJP*JP*JC-((cJCAP*vC*AP*JC)/(v.CFun(t)+vC+cJCAP*AP))-((cJCAC*vC*AC*JC)/(v.CFun(t)+vC+cJCAC*AC))-sC*JC+((100*AC/(50+AC))*wCFun(t))+kCFun(t)
    return(list(c(dAPdt,dACdt,dJPdt,dJCdt)))
  })
}


 ## SI figure 2
#flip scenarios

tstep=1:100

v.PFun=approxfun(x=tstep,y=rep(8,length(tstep)), rule = 2)
v.CFun=approxfun(x=tstep,y=rep(8,length(tstep)), rule = 2)

#### FLIP W/O sp2 harv ####

qEs=seq(0,0.3,length.out = 30)
sto=seq(0,200, length.out = 30)
df2=expand.grid(X=qEs,Y=sto)
df2$AP=numeric(nrow(df2))
df2$AC=numeric(nrow(df2))
df2$JP=numeric(nrow(df2))
df2$JC=numeric(nrow(df2))


for(i in 1:nrow(df2)){
  tstep=1:100
  qEPFun=approxfun(x=tstep,y=rep(df2$X[i], length(tstep)), rule = 2)
  qECFun=approxfun(x=tstep,y=rep(0, length(tstep)), rule = 2)
  kPFun=approxfun(x=tstep,y=rep(df2$Y[i],length(tstep)), rule = 2)
  kCFun=approxfun(x=tstep,y=rep(0,length(tstep)), rule = 2)
  #random % added or subtracted from the deterministic recruitment equation.
  set.seed(1)
  wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
  wC=rtriang(tstep, 0.5, 1.5, 1)
  wPFun=approxfun(x=tstep,y=wP,rule = 2)
  wCFun=approxfun(x=tstep,y=wC,rule = 2)
  
  p=c(sP=0.1,mP=0.1,cJPAP=0.002,cJPAC=0.05,cJPJC=0.003,vP=1,
      sC=0.1,mC=0.1,cJCAC=0.002,cJCAP=0.03,cJCJP=0.003,vC=1)
  y0=c(500,5000,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2w,parms=p)
  df2$AP[i]=sim[nrow(sim)-1,2]
  df2$AC[i]=sim[nrow(sim)-1,3]
  df2$JP[i]=sim[nrow(sim)-1,4]
  df2$JC[i]=sim[nrow(sim)-1,5]
}

df2$diff=df2$AP-df2$AC
df2$ratio=df2$AC/df2$AP
df2$outcome=ifelse(df2$ratio < 0.6,"darkgreen", "darkred")


#### FLIP W/ sp2 harv ####

qEs=seq(0,0.3,length.out = 30)
sto=seq(0,200, length.out = 30)
dfwo2=expand.grid(X=qEs,Y=sto)
dfwo2$AP=numeric(nrow(dfwo2))
dfwo2$AC=numeric(nrow(dfwo2))
dfwo2$JP=numeric(nrow(dfwo2))
dfwo2$JC=numeric(nrow(dfwo2))

for(i in 1:nrow(dfwo2)){
  tstep=1:100
  qEPFun=approxfun(x=tstep,y=rep(dfwo2$X[i], length(tstep)), rule = 2)
  qECFun=approxfun(x=tstep,y=rep(0.1, length(tstep)), rule = 2)
  kPFun=approxfun(x=tstep,y=rep(dfwo2$Y[i],length(tstep)), rule = 2)
  kCFun=approxfun(x=tstep,y=rep(0,length(tstep)), rule = 2)
  #random % added or subtracted from the deterministic recruitment equation.
  set.seed(1)
  wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
  wC=rtriang(tstep, 0.5, 1.5, 1)
  wPFun=approxfun(x=tstep,y=wP,rule = 2)
  wCFun=approxfun(x=tstep,y=wC,rule = 2)
  
  p=c(sP=0.1,mP=0.1,cJPAP=0.002,cJPAC=0.05,cJPJC=0.003,vP=1,
      sC=0.1,mC=0.1,cJCAC=0.002,cJCAP=0.03,cJCJP=0.003,vC=1)
  y0=c(500,5000,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2w,parms=p)
  dfwo2$AP[i]=sim[nrow(sim)-1,2]
  dfwo2$AC[i]=sim[nrow(sim)-1,3]
  dfwo2$JP[i]=sim[nrow(sim)-1,4]
  dfwo2$JC[i]=sim[nrow(sim)-1,5]
}

dfwo2$diff=dfwo2$AP-dfwo2$AC
dfwo2$ratio=dfwo2$AC/dfwo2$AP
dfwo2$outcome=ifelse(dfwo2$ratio < 0.6,"darkgreen", "darkred")
df2$mod=rep("Flip to Preferred, ignore Competitor",nrow(df2))
dfwo2$mod=rep("Flip to Preferred, harv Competitor",nrow(dfwo2))
allSen=rbind(df2,dfwo2)

## use this plot for SI figure 2
si.fig2=ggplot(data=allSen, aes(x=allSen$X,y=allSen$Y,linetype=allSen$mod))+theme_classic()+
  geom_contour(aes(z=allSen$ratio),breaks = c(0.6), color='black', size=1.7)+
  labs(x="Preferred Species Harvest Rate", y="Preferred Species Stocking",linetype="Scenario")+
  theme(legend.position = 'none')+
  geom_segment(aes(x=0.15, y=65, xend=0.2, yend=40), arrow = arrow(type='closed'),size=1)
si.fig2

## SI FIGURE 3
#plot to look at the cost/benefits of stocking or predator reduction for managing a focal species. Flip scenario instead of maintain scenario 

qEs=rep(seq(0,0.3,length.out = 30),3)
sto=rep(seq(0,200, length.out = 30),3)
dfT=expand.grid(X=qEs,Y=sto)
dfT$AP=numeric(nrow(dfT))
dfT$AC=numeric(nrow(dfT))
dfT$JP=numeric(nrow(dfT))
dfT$JC=numeric(nrow(dfT))
dfT$sp1H=c(rep(0.05,30),rep(0.1,30),rep(0.15,30))

for(i in 1:nrow(dfT)){
  tstep=1:100
  qEPFun=approxfun(x=tstep,y=rep(dfT$sp1H[i], length(tstep)))
  qECFun=approxfun(x=tstep,y=rep(dfT$X[i], length(tstep)))
  kPFun=approxfun(x=tstep,y=rep(dfT$Y[i],length(tstep)))
  kCFun=approxfun(x=tstep,y=rep(0,length(tstep)))
  #random % added or subtracted from the deterministic recruitment equation.
  set.seed(1)
  wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
  wC=rtriang(tstep, 0.5, 1.5, 1)
  wPFun=approxfun(x=tstep,y=wP,rule = 2)
  wCFun=approxfun(x=tstep,y=wC,rule = 2)
  
  p=c(sP=0.1,mP=0.1,cJPAP=0.002,cJPAC=0.05,cJPJC=0.003,vP=1,
      sC=0.1,mC=0.1,cJCAC=0.002,cJCAP=0.03,cJCJP=0.003,vC=1)
  y0=c(500,5000,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2w,parms=p)
  dfT$AP[i]=sim[nrow(sim)-1,2]
  dfT$AC[i]=sim[nrow(sim)-1,3]
  dfT$JP[i]=sim[nrow(sim)-1,4]
  dfT$JC[i]=sim[nrow(sim)-1,5]
}

dfT$diff=dfT$AP-dfT$AC
dfT$ratio=dfT$AC/dfT$AP
dfT$qNorm=(dfT$X-min(dfT$X))/(max(dfT$X)-min(dfT$X))
dfT$sNorm=(dfT$Y-min(dfT$Y))/(max(dfT$Y)-min(dfT$Y))
dfT$sp1Norm=((dfT$sp1H-min(dfT$X))/(max(dfT$X)-min(dfT$X))) #putting sp1 harv on same scale as sp2

# SI figure 3
si.fig3=ggplot(data=dfT, aes(x=dfT$X,y=dfT$Y,linetype=as.factor(dfT$sp1H)))+theme_classic()+
  geom_contour(aes(z=dfT$ratio),breaks = c(0.6), color='black', size=1)+
  labs(x="Competitor Species Harvest Rate", y="Preferred Species Stocking",linetype="Preferred Species Harvest")+
  theme(legend.position = 'bottom')
si.fig3
