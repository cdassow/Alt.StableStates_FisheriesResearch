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


## FIGURE 1
#demonstrate alternative stable states
store=data.frame(qEs=seq(0,0.3,length.out=30),A1=0,A2=0,J1=0,J2=0)
y0=c(5000,500,0,0)
tstep=1:100
for(i in 1:nrow(store)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(0.1,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(s1=0.1,m1=0.1,cJ1A1=0.002,cJ1A2=0.05,cJ1J2=0.003,v1=1,
      s2=0.1,m2=0.1,cJ2A2=0.002,cJ2A1=0.03,cJ2J1=0.003,v2=1)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  store$A1[i]=sim[nrow(sim)-1,2]
  store$A2[i]=sim[nrow(sim)-1,3]
  store$J1[i]=sim[nrow(sim)-1,4]
  store$J2[i]=sim[nrow(sim)-1,5]
}
store2=data.frame(qEs=seq(0,0.3,length.out=30),A1=0,A2=0,J1=0,J2=0)
y0=c(500,5000,0,0)
for(i in 1:nrow(store2)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store2$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(0.1,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(s1=0.1,m1=0.1,cJ1A1=0.002,cJ1A2=0.05,cJ1J2=0.003,v1=1,
      s2=0.1,m2=0.1,cJ2A2=0.002,cJ2A1=0.03,cJ2J1=0.003,v2=1)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  store2$A1[i]=sim[nrow(sim)-1,2]
  store2$A2[i]=sim[nrow(sim)-1,3]
  store2$J1[i]=sim[nrow(sim)-1,4]
  store2$J2[i]=sim[nrow(sim)-1,5]
}

## plots
panA=as.data.frame(rbind(cbind(store$qEs,store$A1,rep("A1",nrow(store))),cbind(store$qEs,store$A2,rep("A2",nrow(store)))))
colnames(panA)=c("qEs","Abund","sp")
panA$qEs=as.numeric(panA$qEs);panA$Abund=as.numeric(panA$Abund)
panB=as.data.frame(rbind(cbind(store2$qEs,store2$A1,rep("A1",nrow(store2))),cbind(store2$qEs,store2$A2,rep("A2",nrow(store2)))))
colnames(panB)=c("qEs","Abund","sp")
panB$qEs=as.numeric(panB$qEs);panB$Abund=as.numeric(panB$Abund)

a=ggplot(data = panA,aes(x=qEs,y=Abund,color=sp))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="",labels=c("Preferred Species, initially dominant", "Competitor Species"))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 10),
        legend.position = c(.5,.8))
#a
b=ggplot(data = panB,aes(x=qEs,y=Abund,color=sp))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="",labels=c("Preferred Species", "Competitor Species, initially dominant"))+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 10),
        legend.position = c(.5,.8))
#b
fig1=ggarrange(a,b,ncol=1,labels = c("a","b"))
annotate_figure(fig1,
                left = text_grob("Adult Abundance", rot = 90, size=14),
                bottom = text_grob("Harvest Rate (qE)", size=14))

##FIGURE 2 #####

tstep=1:300

h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)

### MAINTAIN W/O sp2 harv ####
#matrix to hold output, starting with different harvest levels on each species

qEs=seq(0,0.3,length.out = 30)
sto=seq(0,200, length.out = 30)
df=expand.grid(X=qEs,Y=sto)
df$A1=numeric(nrow(df))
df$A2=numeric(nrow(df))
df$J1=numeric(nrow(df))
df$J2=numeric(nrow(df))
minDiff=100

for(i in 1:nrow(df)){
  tstep=1:100
  qE1Fun=approxfun(x=tstep,y=rep(df$X[i], length(tstep)),rule = 2)
  qE2Fun=approxfun(x=tstep,y=rep(0.05, length(tstep)),rule = 2)
  st1Fun=approxfun(x=tstep,y=rep(df$Y[i],length(tstep)),rule = 2)
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  p=c(s1=0.1,m1=0.1,cJ1A1=0.002,cJ1A2=0.05,cJ1J2=0.003,v1=1,
      s2=0.1,m2=0.1,cJ2A2=0.002,cJ2A1=0.03,cJ2J1=0.003,v2=1)
  y0=c(5000,500,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  df$A1[i]=sim[nrow(sim)-1,2]
  df$A2[i]=sim[nrow(sim)-1,3]
  df$J1[i]=sim[nrow(sim)-1,4]
  df$J2[i]=sim[nrow(sim)-1,5]
}

df$diff=df$A1-df$A2
df$ratio=df$A2/df$A1
df$outcome=ifelse(df$ratio < 0.6,"darkgreen", "darkred") # does population fall below 60% of the other



### MAINTAIN W/ sp2 harv ####
#matrix to hold output, starting with different harvest levels on each species

qEs=seq(0,0.3,length.out = 30)
sto=seq(0,200, length.out = 30)
dfwo=expand.grid(X=qEs,Y=sto)
dfwo$A1=numeric(nrow(dfwo))
dfwo$A2=numeric(nrow(dfwo))
dfwo$J1=numeric(nrow(dfwo))
dfwo$J2=numeric(nrow(dfwo))

for(i in 1:nrow(dfwo)){
  tstep=1:100
  qE1Fun=approxfun(x=tstep,y=rep(dfwo$X[i], length(tstep)),rule = 2)
  qE2Fun=approxfun(x=tstep,y=rep(0.15, length(tstep)),rule = 2)
  st1Fun=approxfun(x=tstep,y=rep(dfwo$Y[i],length(tstep)),rule = 2)
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  p=c(s1=0.1,m1=0.1,cJ1A1=0.002,cJ1A2=0.05,cJ1J2=0.003,v1=1,
      s2=0.1,m2=0.1,cJ2A2=0.002,cJ2A1=0.03,cJ2J1=0.003,v2=1)
  y0=c(5000,500,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  dfwo$A1[i]=sim[nrow(sim)-1,2]
  dfwo$A2[i]=sim[nrow(sim)-1,3]
  dfwo$J1[i]=sim[nrow(sim)-1,4]
  dfwo$J2[i]=sim[nrow(sim)-1,5]
}

dfwo$diff=dfwo$A1-dfwo$A2
dfwo$ratio=dfwo$A2/dfwo$A1
dfwo$outcome=ifelse(dfwo$ratio < 0.6,"darkgreen", "darkred") # does population fall below 60% of the other
df$mod=rep("Maintain preferred, ignore competitor",nrow(df))
dfwo$mod=rep("Maintain preferred, harv competitor",nrow(dfwo))
allSen=rbind(df,dfwo)
allSen$percDiff=(abs(allSen$A1-allSen$A2)/((allSen$A1+allSen$A2)/2))*100


# vzI=ggplot(data=allSen, aes(x=allSen$X,y=allSen$Y,linetype=allSen$mod))+
#   geom_contour(aes(z=allSen$diff), breaks=c(100), color='black', size=1.7)+
#   labs(x="Preferred Species Harvest Rate", y="Preferred Stocking",linetype="Scenario")+
#   theme(legend.position = 'right') + theme(legend.position = 'none')
# vzI
# 
# 
# ggplot(allSen)+theme_classic()+
#   geom_tile(aes(x=X, y=Y, fill=outcome))+
#   facet_wrap(~mod)+
#   coord_cartesian(ylim=c(0,200))

### use this one
ggplot(data=allSen, aes(x=allSen$X,y=allSen$Y,linetype=allSen$mod))+theme_classic()+
  geom_contour(aes(z=allSen$ratio),breaks = c(0.6), color='black', size=1.7)+
  geom_segment(aes(x=0.20, y=30, xend=0.25, yend=10), arrow = arrow(type='closed'),size=1)+
  labs(x="Preferred Species Harvest Rate", y="Preferred Species Stocking",linetype="Scenario")+
  theme(legend.position = 'right') + theme(legend.position = 'none')


### FIGURE 3 ####
#figure 3, isoclines at different harvests

#plot to look at the cost/benefits of stocking or predator reduction for managing a focal species. 

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
  y0=c(5000,500,0,0)
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

vzT=ggplot(data=dfT, aes(x=dfT$X,y=dfT$Y,linetype=as.factor(dfT$sp1H)))+theme_classic()+
  geom_contour(aes(z=dfT$ratio),breaks = c(0.6), color='black', size=1)+
  labs(x="Competitor Species Harvest Rate", y="Preferred Species Stocking",linetype="Preferred Species Harvest")+
  theme(legend.position = c(.75,.75))+
  scale_linetype_manual(values = c("solid","dashed","twodash"))
vzT


### FIGURE 4, DELAY A TRANSITION ####
#same function as used in other figs but with option to change fecundity now
simBiggsR.a<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
    dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
    dJ1dt=-cJ1J2*J2*J1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))-((cJ1A1*v1*A1*J1)/(h1Fun(t)+v1+cJ1A1*A1))-s1*J1+((100*dFun(t))*A1/(50+A1))+st1Fun(t)
    dJ2dt=-cJ2J1*J1*J2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))-((cJ2A2*v2*A2*J2)/(h2Fun(t)+v2+cJ2A2*A2))-s2*J2+(100*A2/(50+A2))+st2Fun(t)
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}

#running to sp1 (prefered species) dominating before habitat decline
tstep=1:100
h1Fun=approxfun(x=tstep,y=c(rep(8,length(tstep))),rule = 2)
h2Fun=approxfun(x=tstep,y=c(rep(8,length(tstep))),rule = 2)
qE1Fun=approxfun(x=tstep,y=c(rep(0,length(tstep))),rule = 2)
qE2Fun=approxfun(x=tstep,y=c(rep(0,length(tstep))),rule = 2)
st1Fun=approxfun(x=tstep,y=c(rep(0,length(tstep))),rule = 2)
st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
dFun=approxfun(x=tstep,y=rep(1,length(tstep)),rule = 2)

p=c(s1=0.1,m1=0.1,cJ1A1=0.001,cJ1A2=0.05,cJ1J2=0.003,v1=1,
    s2=0.1,m2=0.1,cJ2A2=0.001,cJ2A1=0.03,cJ2J1=0.003,v2=1)
y0=c(5000,500,0,0)
simPre=ode(y=y0,times=tstep,func=simBiggsR.a,parms=p)

#Fecundity decline -> what will happen if no action taken panel a
tstep=1:100
h1Fun=approxfun(x=tstep,y=c(rep(8,100)),rule = 2)
h2Fun=approxfun(x=tstep,y=c(rep(8,100)),rule = 2)
qE1Fun=approxfun(x=tstep,y=c(rep(0,100)),rule = 2)
qE2Fun=approxfun(x=tstep,y=c(rep(0,100)),rule = 2)
st1Fun=approxfun(x=tstep,y=c(rep(0,100)),rule = 2)
st2Fun=approxfun(x=tstep,y=c(rep(0,100)),rule = 2)
dFun=approxfun(x=tstep,y=c(seq(1,.01,length.out = 20),rep(.01,80)),rule = 2)

p=c(s1=0.1,m1=0.1,cJ1A1=0.002,cJ1A2=0.05,cJ1J2=0.003,v1=1,
    s2=0.1,m2=0.1,cJ2A2=0.002,cJ2A1=0.03,cJ2J1=0.003,v2=1)
#y0=simPre[nrow(simPre)-1,2:5]
sim=ode(y=y0,times=tstep,func=simBiggsR.a,parms=p)

#stocking delays the transition to sp2 (competitor), panel c
tstep=1:100
h1Fun=approxfun(x=tstep,y=c(rep(8,100)),rule = 2)
h2Fun=approxfun(x=tstep,y=c(rep(8,100)),rule = 2)
qE1Fun=approxfun(x=tstep,y=c(rep(0,100)),rule = 2)
qE2Fun=approxfun(x=tstep,y=c(rep(0,100)),rule = 2)
st1Fun=approxfun(x=tstep,y=c(rep(50,100)),rule = 2)
st2Fun=approxfun(x=tstep,y=c(rep(0,100)),rule = 2) 
dFun=approxfun(x=tstep,y=c(seq(1,0.01,length.out = 20),rep(.01,80)),rule = 2)

simS=ode(y=y0,times=tstep,func=simBiggsR.a,parms=p)

#harvesting delays the transition to sp2 - no stocking, panel b
tstep=1:100
h1Fun=approxfun(x=tstep,y=c(rep(8,100)),rule = 2)
h2Fun=approxfun(x=tstep,y=c(rep(8,100)),rule = 2)
qE1Fun=approxfun(x=tstep,y=c(rep(0,100)),rule = 2)
qE2Fun=approxfun(x=tstep,y=c(rep(0.15,100)),rule = 2)
st1Fun=approxfun(x=tstep,y=c(rep(0,100)),rule = 2)
st2Fun=approxfun(x=tstep,y=c(rep(0,100)),rule = 2) 
dFun=approxfun(x=tstep,y=c(seq(1,0.01,length.out = 20),rep(.01,80)),rule = 2)

simH=ode(y=y0,times=tstep,func=simBiggsR.a,parms=p)

#harvesting and stocking delays the transition to sp2 indefinitely, panel d
tstep=1:100
h1Fun=approxfun(x=tstep,y=c(rep(8,100)),rule = 2)
h2Fun=approxfun(x=tstep,y=c(rep(8,100)),rule = 2)
qE1Fun=approxfun(x=tstep,y=c(rep(0,100)),rule = 2)
qE2Fun=approxfun(x=tstep,y=c(rep(0.15,100)),rule = 2)
st1Fun=approxfun(x=tstep,y=c(rep(50,100)),rule = 2)
st2Fun=approxfun(x=tstep,y=c(rep(0,100)),rule = 2) 
dFun=approxfun(x=tstep,y=c(seq(1,0.01,length.out = 20),rep(.01,80)),rule = 2)

simB=ode(y=y0,times=tstep,func=simBiggsR.a,parms=p)

#plotting

#do nothing
panA=as.data.frame(rbind(cbind(sim[,c(1,2)],rep("A1",nrow(sim))),cbind(sim[,c(1,3)],rep("A2",nrow(sim)))))
colnames(panA)=c("Time","Abund","sp")
panA$Abund=as.numeric(panA$Abund);panA$Time=as.numeric(panA$Time)
transTime.a=which(abs(sim[,2]-sim[,3])==min(abs(sim[,2]-sim[,3])))
transN.a=sim[transTime.a,2]

#harvest sp2
panB=as.data.frame(rbind(cbind(simH[,c(1,2)],rep("A1",nrow(simH))),cbind(simH[,c(1,3)],rep("A2",nrow(simH)))))
colnames(panB)=c("Time","Abund","sp")
panB$Time=as.numeric(panB$Time);panB$Abund=as.numeric(panB$Abund)
transTime.b=which(abs(simH[,2]-simH[,3])==min(abs(simH[,2]-simH[,3])))
transN.b=sim[transTime.b,2]

#stock sp1
panC=as.data.frame(rbind(cbind(simS[,c(1,2)],rep("A1",nrow(simS))),cbind(simS[,c(1,3)],rep("A2",nrow(simS)))))
colnames(panC)=c("Time","Abund","sp")
panC$Abund=as.numeric(panC$Abund);panC$Time=as.numeric(panC$Time)
transTime.c=which(abs(simS[,2]-simS[,3])==min(abs(simS[,2]-simS[,3])))
transN.c=sim[transTime.c,2]

#stock and harvest
panD=as.data.frame(rbind(cbind(simB[,c(1,2)],rep("A1",nrow(simB))),cbind(simB[,c(1,3)],rep("A2",nrow(simB)))))
colnames(panD)=c("Time","Abund","sp")
panD$Abund=as.numeric(panD$Abund);panD$Time=as.numeric(panD$Time)
transTime.d=which(abs(simB[,2]-simB[,3])==min(abs(simB[,2]-simB[,3])))
transN.d=sim[transTime.d,2]

a4=ggplot(data = panA,aes(x=Time,y=Abund,color=sp))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())+
 geom_vline(xintercept = transTime.a)+
  ylim(0,2000)
b4=ggplot(data = panB,aes(x=Time,y=Abund,color=sp))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  geom_vline(xintercept = transTime.b)+
  ylim(-1,2000)
c4=ggplot(data = panC,aes(x=Time,y=Abund,color=sp))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  geom_vline(xintercept = transTime.c)+
  ylim(-1,2000)
d4=ggplot(data = panD,aes(x=Time,y=Abund,color=sp))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  ylim(-1,2000)

fig4=ggarrange(a4,b4,c4,d4,labels = c("a","b","c","d"),common.legend = T,legend = "top",hjust = c(0.5,1,0.5,1))
annotate_figure(fig4,
                left = text_grob("Adult Abundance", rot = 90, vjust=0),
                bottom = text_grob("Time"))
