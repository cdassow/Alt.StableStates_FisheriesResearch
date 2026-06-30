# script to repeat analysis but with stochasticity added to recruitment in the model.
## CJD 6.11.26

rm(list=ls())
library(deSolve)
library(ggplot2)
library(ggpubr)
library(extraDistr)
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


## FIGURE 1
#demonstrate alternative stable states
store=data.frame(qEs=seq(0,0.3,length.out=30),AP=0,AC=0,JP=0,JC=0)
y0=c(5000,500,0,0) # c(AP, AC, JP, JC)
tstep=1:100
for(i in 1:nrow(store)){
  qEPFun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))), rule = 2)
  qECFun=approxfun(x=tstep,y=rep(0.1,length(tstep)), rule = 2)
  v.PFun=approxfun(x=tstep,y=rep(8,length(tstep)), rule = 2)
  v.CFun=approxfun(x=tstep,y=rep(8,length(tstep)), rule = 2)
  kPFun=approxfun(x=tstep,y=rep(0,length(tstep)), rule = 2)
  kCFun=approxfun(x=tstep,y=rep(0,length(tstep)), rule = 2)

  p=c(sP=0.1,mP=0.1,cJPAP=0.002,cJPAC=0.05,cJPJC=0.003,vP=1,
      sC=0.1,mC=0.1,cJCAC=0.002,cJCAP=0.03,cJCJP=0.003,vC=1)
  sim=data.frame(time=NA, AP=NA, AC=NA, JP=NA, JC=NA,modelRun=NA) # data frame to store each of the 10 runs in
  set.seed(1)
  for(r in 1:10){
    #random % added or subtracted from the deterministic recruitment equation.
    wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
    wC=rtriang(tstep, 0.5, 1.5, 1)
    wPFun=approxfun(x=tstep,y=wP,rule = 2)
    wCFun=approxfun(x=tstep,y=wC,rule = 2)
    
    tsim=ode(y=y0,times=tstep,func=simBiggsQ2w,parms=p)
    tsim=as.data.frame(cbind(tsim, rep(r,length(tstep))))
    colnames(tsim)=colnames(sim)
    sim=rbind(sim,tsim)
  }
  sim=sim[!is.na(sim$time),] # removing row of NAs used to initialize dataframe
  store$AP[i]=mean(sim$AP[sim$time==length(tstep)-1])# taking the mean of all 10 runs
  store$AC[i]=mean(sim$AC[sim$time==length(tstep)-1])
  store$JP[i]=mean(sim$JP[sim$time==length(tstep)-1])
  store$JC[i]=mean(sim$JC[sim$time==length(tstep)-1])
}
store2=data.frame(qEs=seq(0,0.3,length.out=30),AP=0,AC=0,JP=0,JC=0)
y0=c(500,5000,0,0)
for(i in 1:nrow(store2)){
  qEPFun=approxfun(x=tstep,y=c(rep(store2$qEs[i],length(tstep))), rule = 2)
  qECFun=approxfun(x=tstep,y=rep(0.1,length(tstep)), rule = 2)
  v.PFun=approxfun(x=tstep,y=rep(8,length(tstep)), rule = 2)
  v.CFun=approxfun(x=tstep,y=rep(8,length(tstep)), rule = 2)
  kPFun=approxfun(x=tstep,y=rep(0,length(tstep)), rule = 2)
  kCFun=approxfun(x=tstep,y=rep(0,length(tstep)), rule = 2)

  p=c(sP=0.1,mP=0.1,cJPAP=0.002,cJPAC=0.05,cJPJC=0.003,vP=1,
      sC=0.1,mC=0.1,cJCAC=0.002,cJCAP=0.03,cJCJP=0.003,vC=1)
  sim=data.frame(time=NA, AP=NA, AC=NA, JP=NA, JC=NA,modelRun=NA) # data frame to store each of the 10 runs in
  set.seed(1)
  for(r in 1:10){
    #random % added or subtracted from the deterministic recruitment equation.
    wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
    wC=rtriang(tstep, 0.5, 1.5, 1)
    wPFun=approxfun(x=tstep,y=wP,rule = 2)
    wCFun=approxfun(x=tstep,y=wC,rule = 2)
    
    tsim=ode(y=y0,times=tstep,func=simBiggsQ2w,parms=p)
    tsim=as.data.frame(cbind(tsim, rep(r,length(tstep))))
    colnames(tsim)=colnames(sim)
    sim=rbind(sim,tsim)
  }
  sim=sim[!is.na(sim$time),] # removing row of NAs used to initialize dataframe
  store2$AP[i]=mean(sim$AP[sim$time==length(tstep)-1])# taking the mean of all 10 runs
  store2$AC[i]=mean(sim$AC[sim$time==length(tstep)-1])
  store2$JP[i]=mean(sim$JP[sim$time==length(tstep)-1])
  store2$JC[i]=mean(sim$JC[sim$time==length(tstep)-1])
}

## plots
panA.w=as.data.frame(rbind(cbind(store$qEs,store$AP,rep("AP",nrow(store))),cbind(store$qEs,store$AC,rep("AC",nrow(store)))))
colnames(panA.w)=c("qEs","Abund","sp")
panA.w$qEs=as.numeric(panA.w$qEs);panA.w$Abund=as.numeric(panA.w$Abund)
panB.w=as.data.frame(rbind(cbind(store2$qEs,store2$AP,rep("AP",nrow(store2))),cbind(store2$qEs,store2$AC,rep("AC",nrow(store2)))))
colnames(panB.w)=c("qEs","Abund","sp")
panB.w$qEs=as.numeric(panB.w$qEs);panB.w$Abund=as.numeric(panB.w$Abund)

a=ggplot(data = panA.w,aes(x=qEs,y=Abund,color=sp))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="",labels=c("Preferred Species, initially dominant", "Competitor Species"))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 10),
        legend.position = c(.5,.85))
#a
b=ggplot(data = panB.w,aes(x=qEs,y=Abund,color=sp))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="",labels=c("Preferred Species", "Competitor Species, initially dominant"))+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 10),
        legend.position = c(.5,.85))
#b
alt.fig1=ggarrange(a,b,ncol=1,labels = c("a","b"))
alt.fig1 = annotate_figure(alt.fig1,
                left = text_grob("Adult Abundance", rot = 90, size=14),
                bottom = text_grob("Harvest Rate (qE)", size=14))
alt.fig1
ggsave('alt.fig1_revision.jpeg', plot = alt.fig1, dpi = 700)


##FIGURE 2 #####

tstep=1:100

v.PFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
v.CFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)

### MAINTAIN W/O sp2 harv ####
#matrix to hold output, starting with different harvest levels on each species

qEs=seq(0,0.3,length.out = 30)
sto=seq(0,200, length.out = 30)
df=expand.grid(X=qEs,Y=sto)
df$AP=numeric(nrow(df))
df$AC=numeric(nrow(df))
df$JP=numeric(nrow(df))
df$JC=numeric(nrow(df))

for(i in 1:nrow(df)){
  tstep=1:100
  qEPFun=approxfun(x=tstep,y=rep(df$X[i], length(tstep)),rule = 2)
  qECFun=approxfun(x=tstep,y=rep(0.05, length(tstep)),rule = 2)
  kPFun=approxfun(x=tstep,y=rep(df$Y[i],length(tstep)),rule = 2)
  kCFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)

  p=c(sP=0.1,mP=0.1,cJPAP=0.002,cJPAC=0.05,cJPJC=0.003,vP=1,
      sC=0.1,mC=0.1,cJCAC=0.002,cJCAP=0.03,cJCJP=0.003,vC=1)
  y0=c(5000,500,0,0)
  sim=data.frame(time=NA, AP=NA, AC=NA, JP=NA, JC=NA,modelRun=NA) # data frame to store each of the 10 runs in
  set.seed(1)
  for(r in 1:10){
    #random % added or subtracted from the deterministic recruitment equation.
    wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
    wC=rtriang(tstep, 0.5, 1.5, 1)
    wPFun=approxfun(x=tstep,y=wP,rule = 2)
    wCFun=approxfun(x=tstep,y=wC,rule = 2)
    
    tsim=ode(y=y0,times=tstep,func=simBiggsQ2w,parms=p)
    tsim=as.data.frame(cbind(tsim, rep(r,length(tstep))))
    colnames(tsim)=colnames(sim)
    sim=rbind(sim,tsim)
  }
  sim=sim[!is.na(sim$time),] # removing row of NAs used to initialize dataframe
  df$AP[i]=mean(sim$AP[sim$time==length(tstep)-1])# taking the mean of all 10 runs
  df$AC[i]=mean(sim$AC[sim$time==length(tstep)-1])
  df$JP[i]=mean(sim$JP[sim$time==length(tstep)-1])
  df$JC[i]=mean(sim$JC[sim$time==length(tstep)-1])
}

df$diff=df$AP-df$AC
df$ratio=df$AC/df$AP
df$outcome=ifelse(df$ratio <= 0.6,"good", "bad") # does population fall below 60% of the other



### MAINTAIN W/ sp2 harv ####
#matrix to hold output, starting with different harvest levels on each species

qEs=seq(0,0.3,length.out = 30)
sto=seq(0,200, length.out = 30)
dfwo=expand.grid(X=qEs,Y=sto)
dfwo$AP=numeric(nrow(dfwo))
dfwo$AC=numeric(nrow(dfwo))
dfwo$JP=numeric(nrow(dfwo))
dfwo$JC=numeric(nrow(dfwo))

for(i in 1:nrow(dfwo)){
  tstep=1:100
  qEPFun=approxfun(x=tstep,y=rep(dfwo$X[i], length(tstep)),rule = 2)
  qECFun=approxfun(x=tstep,y=rep(0.15, length(tstep)),rule = 2)
  kPFun=approxfun(x=tstep,y=rep(dfwo$Y[i],length(tstep)),rule = 2)
  kCFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)

  p=c(sP=0.1,mP=0.1,cJPAP=0.002,cJPAC=0.05,cJPJC=0.003,vP=1,
      sC=0.1,mC=0.1,cJCAC=0.002,cJCAP=0.03,cJCJP=0.003,vC=1)
  y0=c(5000,500,0,0)
  sim=data.frame(time=NA, AP=NA, AC=NA, JP=NA, JC=NA,modelRun=NA) # data frame to store each of the 10 runs in
  set.seed(1)
  for(r in 1:10){
    #random % added or subtracted from the deterministic recruitment equation.
    wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
    wC=rtriang(tstep, 0.5, 1.5, 1)
    wPFun=approxfun(x=tstep,y=wP,rule = 2)
    wCFun=approxfun(x=tstep,y=wC,rule = 2)
    
    tsim=ode(y=y0,times=tstep,func=simBiggsQ2w,parms=p)
    tsim=as.data.frame(cbind(tsim, rep(r,length(tstep))))
    colnames(tsim)=colnames(sim)
    sim=rbind(sim,tsim)
  }
  sim=sim[!is.na(sim$time),] # removing row of NAs used to initialize dataframe
  dfwo$AP[i]=mean(sim$AP[sim$time==length(tstep)-1])# taking the mean of all 10 runs
  dfwo$AC[i]=mean(sim$AC[sim$time==length(tstep)-1])
  dfwo$JP[i]=mean(sim$JP[sim$time==length(tstep)-1])
  dfwo$JC[i]=mean(sim$JC[sim$time==length(tstep)-1])
}

dfwo$diff=dfwo$AP-dfwo$AC
dfwo$ratio=dfwo$AC/dfwo$AP
dfwo$outcome=ifelse(dfwo$ratio <= 0.6,"good", "bad") # does population fall below 60% of the other
df$mod=rep("Maintain preferred, ignore competitor",nrow(df))
dfwo$mod=rep("Maintain preferred, harv competitor",nrow(dfwo))
allSen.w=rbind(df,dfwo)
allSen.w$percDiff=(abs(allSen.w$AP-allSen.w$AC)/((allSen.w$AP+allSen.w$AC)/2))*100


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
alt.fig2=ggplot(data=allSen.w, aes(x=allSen.w$X,y=allSen.w$Y,linetype=allSen.w$mod))+theme_classic()+
  geom_contour(aes(z=allSen.w$ratio),breaks = c(0.6), color='black', size=1.7)+
  geom_segment(aes(x=0.20, y=30, xend=0.25, yend=10), arrow = arrow(type='closed'),size=1)+
  labs(x="Preferred Species Harvest Rate", y="Preferred Species Stocking",linetype="Scenario")+
  theme(legend.position = 'none')+coord_cartesian(ylim=c(0,100), xlim=c(0,0.3))
alt.fig2
ggsave('alt.fig2_revision.jpeg', plot = alt.fig2, dpi=700)


### FIGURE 3 ####
#figure 3, isoclines at different harvests

#plot to look at the cost/benefits of stocking or predator reduction for managing a focal species. 

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
  qEPFun=approxfun(x=tstep,y=rep(dfT$sp1H[i], length(tstep)), rule = 2)
  qECFun=approxfun(x=tstep,y=rep(dfT$X[i], length(tstep)), rule = 2)
  kPFun=approxfun(x=tstep,y=rep(dfT$Y[i],length(tstep)), rule = 2)
  kCFun=approxfun(x=tstep,y=rep(0,length(tstep)), rule = 2)

  p=c(sP=0.1,mP=0.1,cJPAP=0.002,cJPAC=0.05,cJPJC=0.003,vP=1,
      sC=0.1,mC=0.1,cJCAC=0.002,cJCAP=0.03,cJCJP=0.003,vC=1)
  y0=c(5000,500,0,0)
  sim=data.frame(time=NA, AP=NA, AC=NA, JP=NA, JC=NA,modelRun=NA) # data frame to store each of the 10 runs in
  set.seed(1)
  for(r in 1:10){
    #random % added or subtracted from the deterministic recruitment equation.
    wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
    wC=rtriang(tstep, 0.5, 1.5, 1)
    wPFun=approxfun(x=tstep,y=wP,rule = 2)
    wCFun=approxfun(x=tstep,y=wC,rule = 2)
    
    tsim=ode(y=y0,times=tstep,func=simBiggsQ2w,parms=p)
    tsim=as.data.frame(cbind(tsim, rep(r,length(tstep))))
    colnames(tsim)=colnames(sim)
    sim=rbind(sim,tsim)
  }
  sim=sim[!is.na(sim$time),] # removing row of NAs used to initialize dataframe
  dfT$AP[i]=mean(sim$AP[sim$time==length(tstep)-1])# taking the mean of all 10 runs
  dfT$AC[i]=mean(sim$AC[sim$time==length(tstep)-1])
  dfT$JP[i]=mean(sim$JP[sim$time==length(tstep)-1])
  dfT$JC[i]=mean(sim$JC[sim$time==length(tstep)-1])
}

dfT$diff=dfT$AP-dfT$AC
dfT$ratio=dfT$AC/dfT$AP

alt.fig3=ggplot(data=dfT, aes(x=dfT$X,y=dfT$Y,linetype=as.factor(dfT$sp1H)))+theme_classic()+
  geom_contour(aes(z=dfT$ratio),breaks = c(0.6), color='black', size=1)+
  labs(x="Competitor Species Harvest Rate", y="Preferred Species Stocking",linetype="Preferred Species Harvest")+
  theme(legend.position = c(.75,.75), legend.key.width = unit(3, 'line'))+
  scale_linetype_manual(values = c("solid","dashed","twodash"))+
  coord_cartesian(xlim=c(0,0.15), ylim=c(0,80))
alt.fig3
ggsave('alt.fig3_revision.jpeg', plot = alt.fig3, dpi=700)

### FIGURE 4, DELAY A TRANSITION ####
#same function as used in other figs but with option to change fecundity now
simBiggsR.aw<-function(t,y,params){
  AP<-y[1]
  AC<-y[2]
  JP<-y[3]
  JC<-y[4]
  with(as.list(params),{
    dAPdt=-qEPFun(t)*AP-mP*AP+sP*JP
    dACdt=-qECFun(t)*AC-mC*AC+sC*JC
    dJPdt=-cJPJC*JC*JP-((cJPAC*vP*AC*JP)/(v.PFun(t)+vP+cJPAC*AC))-((cJPAP*vP*AP*JP)/(v.PFun(t)+vP+cJPAP*AP))-sP*JP+(((100*dFun(t))*AP/(50+AP))*wPFun(t))+kPFun(t)
    dJCdt=-cJCJP*JP*JC-((cJCAP*vC*AP*JC)/(v.CFun(t)+vC+cJCAP*AP))-((cJCAC*vC*AC*JC)/(v.CFun(t)+vC+cJCAC*AC))-sC*JC+((100*AC/(50+AC))*wCFun(t))+kCFun(t)
    return(list(c(dAPdt,dACdt,dJPdt,dJCdt)))
  })
}

#running to prefered species dominating before habitat decline
# not used in final analysis
# tstep=1:100
# v.PFun=approxfun(x=tstep,y=c(rep(8,length(tstep))),rule = 2)
# v.CFun=approxfun(x=tstep,y=c(rep(8,length(tstep))),rule = 2)
# qEPFun=approxfun(x=tstep,y=c(rep(0,length(tstep))),rule = 2)
# qECFun=approxfun(x=tstep,y=c(rep(0,length(tstep))),rule = 2)
# kPFun=approxfun(x=tstep,y=c(rep(0,length(tstep))),rule = 2)
# kCFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
# dFun=approxfun(x=tstep,y=rep(1,length(tstep)),rule = 2)
# 
# p=c(sP=0.1,mP=0.1,cJPAP=0.001,cJPAC=0.05,cJPJC=0.003,vP=1,
#     sC=0.1,mC=0.1,cJCAC=0.001,cJCAP=0.03,cJCJP=0.003,vC=1)
# y0=c(5000,500,0,0)
# simPre=ode(y=y0,times=tstep,func=simBiggsR.a,parms=p)

#Fecundity decline -> what will happen if no action taken panel a
tstep=1:100
v.PFun=approxfun(x=tstep,y=c(rep(8,100)),rule = 2)
v.CFun=approxfun(x=tstep,y=c(rep(8,100)),rule = 2)
qEPFun=approxfun(x=tstep,y=c(rep(0,100)),rule = 2)
qECFun=approxfun(x=tstep,y=c(rep(0,100)),rule = 2)
kPFun=approxfun(x=tstep,y=c(rep(0,100)),rule = 2)
kCFun=approxfun(x=tstep,y=c(rep(0,100)),rule = 2)
dFun=approxfun(x=tstep,y=c(seq(1,.01,length.out = 20),rep(.01,80)),rule = 2)

p=c(sP=0.1,mP=0.1,cJPAP=0.002,cJPAC=0.05,cJPJC=0.003,vP=1,
    sC=0.1,mC=0.1,cJCAC=0.002,cJCAP=0.03,cJCJP=0.003,vC=1)
#y0=simPre[nrow(simPre)-1,2:5]
sim=data.frame(time=NA, AP=NA, AC=NA, JP=NA, JC=NA,modelRun=NA) # data frame to store each of the 10 runs in
set.seed(1)
for(r in 1:10){
  #random % added or subtracted from the deterministic recruitment equation.
  wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
  wC=rtriang(tstep, 0.5, 1.5, 1)
  wPFun=approxfun(x=tstep,y=wP,rule = 2)
  wCFun=approxfun(x=tstep,y=wC,rule = 2)
  
  tsim=ode(y=y0,times=tstep,func=simBiggsR.aw,parms=p)
  tsim=as.data.frame(cbind(tsim, rep(r,length(tstep))))
  colnames(tsim)=colnames(sim)
  sim=rbind(sim,tsim)
}
sim=sim[!is.na(sim$time),] # removing row of NAs used to initialize dataframe
sim.a=data.frame(time=tstep,
                 AP=NA,
                 AC=NA,
                 JP=NA,
                 JC=NA)
for(i in 1:nrow(sim.a)){
  sim.a$AP[i]=mean(sim$AP[sim$time==sim.a$time[i]])
  sim.a$AC[i]=mean(sim$AC[sim$time==sim.a$time[i]])
  sim.a$JP[i]=mean(sim$JP[sim$time==sim.a$time[i]])
  sim.a$JC[i]=mean(sim$JC[sim$time==sim.a$time[i]])
  
}

#stocking delays the transition to competitor, panel c
tstep=1:100
v.PFun=approxfun(x=tstep,y=c(rep(8,100)),rule = 2)
v.CFun=approxfun(x=tstep,y=c(rep(8,100)),rule = 2)
qEPFun=approxfun(x=tstep,y=c(rep(0,100)),rule = 2)
qECFun=approxfun(x=tstep,y=c(rep(0,100)),rule = 2)
kPFun=approxfun(x=tstep,y=c(rep(50,100)),rule = 2) # stocking the preferred species
kCFun=approxfun(x=tstep,y=c(rep(0,100)),rule = 2) 
dFun=approxfun(x=tstep,y=c(seq(1,0.01,length.out = 20),rep(.01,80)),rule = 2)

simS=data.frame(time=NA, AP=NA, AC=NA, JP=NA, JC=NA,modelRun=NA) # data frame to store each of the 10 runs in
set.seed(1)
for(r in 1:10){
  #random % added or subtracted from the deterministic recruitment equation.
  wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
  wC=rtriang(tstep, 0.5, 1.5, 1)
  wPFun=approxfun(x=tstep,y=wP,rule = 2)
  wCFun=approxfun(x=tstep,y=wC,rule = 2)
  
  tsim=ode(y=y0,times=tstep,func=simBiggsR.aw,parms=p)
  tsim=as.data.frame(cbind(tsim, rep(r,length(tstep))))
  colnames(tsim)=colnames(simS)
  simS=rbind(simS,tsim)
}
simS=simS[!is.na(simS$time),] # removing row of NAs used to initialize dataframe

sim.c=data.frame(time=tstep,
                 AP=NA,
                 AC=NA,
                 JP=NA,
                 JC=NA)
for(i in 1:nrow(sim.c)){
  sim.c$AP[i]=mean(simS$AP[simS$time==sim.c$time[i]])
  sim.c$AC[i]=mean(simS$AC[simS$time==sim.c$time[i]])
  sim.c$JP[i]=mean(simS$JP[simS$time==sim.c$time[i]])
  sim.c$JC[i]=mean(simS$JC[simS$time==sim.c$time[i]])
  
}

#harvesting delays the transition to competitor - no stocking, panel b
tstep=1:100
v.PFun=approxfun(x=tstep,y=c(rep(8,100)),rule = 2)
v.CFun=approxfun(x=tstep,y=c(rep(8,100)),rule = 2)
qEPFun=approxfun(x=tstep,y=c(rep(0,100)),rule = 2)
qECFun=approxfun(x=tstep,y=c(rep(0.15,100)),rule = 2) # harvesting the competitor species
kPFun=approxfun(x=tstep,y=c(rep(0,100)),rule = 2)
kCFun=approxfun(x=tstep,y=c(rep(0,100)),rule = 2) 
dFun=approxfun(x=tstep,y=c(seq(1,0.01,length.out = 20),rep(.01,80)),rule = 2)

simH=data.frame(time=NA, AP=NA, AC=NA, JP=NA, JC=NA,modelRun=NA) # data frame to store each of the 10 runs in
set.seed(1)
for(r in 1:10){
  #random % added or subtracted from the deterministic recruitment equation.
  wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
  wC=rtriang(tstep, 0.5, 1.5, 1)
  wPFun=approxfun(x=tstep,y=wP,rule = 2)
  wCFun=approxfun(x=tstep,y=wC,rule = 2)
  
  tsim=ode(y=y0,times=tstep,func=simBiggsR.aw,parms=p)
  tsim=as.data.frame(cbind(tsim, rep(r,length(tstep))))
  colnames(tsim)=colnames(simH)
  simH=rbind(simH,tsim)
}
simH=simH[!is.na(simH$time),] # removing row of NAs used to initialize dataframe

sim.b=data.frame(time=tstep,
                 AP=NA,
                 AC=NA,
                 JP=NA,
                 JC=NA)
for(i in 1:nrow(sim.b)){
  sim.b$AP[i]=mean(simH$AP[simH$time==sim.b$time[i]])
  sim.b$AC[i]=mean(simH$AC[simH$time==sim.b$time[i]])
  sim.b$JP[i]=mean(simH$JP[simH$time==sim.b$time[i]])
  sim.b$JC[i]=mean(simH$JC[simH$time==sim.b$time[i]])
  
}
#harvesting and stocking delays the transition to competitor indefinitely, panel d
tstep=1:100
v.PFun=approxfun(x=tstep,y=c(rep(8,100)),rule = 2)
v.CFun=approxfun(x=tstep,y=c(rep(8,100)),rule = 2)
qEPFun=approxfun(x=tstep,y=c(rep(0,100)),rule = 2)
qECFun=approxfun(x=tstep,y=c(rep(0.15,100)),rule = 2) # harvesting the competitor
kPFun=approxfun(x=tstep,y=c(rep(50,100)),rule = 2) # stocking the preferred
kCFun=approxfun(x=tstep,y=c(rep(0,100)),rule = 2) 
dFun=approxfun(x=tstep,y=c(seq(1,0.01,length.out = 20),rep(.01,80)),rule = 2)

simB=data.frame(time=NA, AP=NA, AC=NA, JP=NA, JC=NA,modelRun=NA) # data frame to store each of the 10 runs in
set.seed(1)
for(r in 1:10){
  #random % added or subtracted from the deterministic recruitment equation.
  wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
  wC=rtriang(tstep, 0.5, 1.5, 1)
  wPFun=approxfun(x=tstep,y=wP,rule = 2)
  wCFun=approxfun(x=tstep,y=wC,rule = 2)
  
  tsim=ode(y=y0,times=tstep,func=simBiggsR.aw,parms=p)
  tsim=as.data.frame(cbind(tsim, rep(r,length(tstep))))
  colnames(tsim)=colnames(simB)
  simB=rbind(simB,tsim)
}
simB=simB[!is.na(simB$time),] # removing row of NAs used to initialize dataframe

sim.d=data.frame(time=tstep,
                 AP=NA,
                 AC=NA,
                 JP=NA,
                 JC=NA)
for(i in 1:nrow(sim.d)){
  sim.d$AP[i]=mean(simB$AP[simB$time==sim.d$time[i]])
  sim.d$AC[i]=mean(simB$AC[simB$time==sim.d$time[i]])
  sim.d$JP[i]=mean(simB$JP[simB$time==sim.d$time[i]])
  sim.d$JC[i]=mean(simB$JC[simB$time==sim.d$time[i]])
  
}

#plotting

#do nothing
tp=cbind(sim.a[,c(1,2)],rep('Preferred Species',nrow(sim.a)),deparse.level = 0)
colnames(tp)=c("Time","Abund","sp")
tc=cbind(sim.a[,c(1,3)],rep('Competitor Species',nrow(sim.a)),deparse.level = 0)
colnames(tc)=c("Time","Abund","sp")
panA=as.data.frame(rbind(tp,tc))

transTime.a=which(abs(sim.a[,2]-sim.a[,3])==min(abs(sim.a[,2]-sim.a[,3])))
transN.a=sim.a[transTime.a,2]

#harvest sp2
tp=cbind(sim.b[,c(1,2)],rep('Preferred Species',nrow(sim.b)),deparse.level = 0)
colnames(tp)=c("Time","Abund","sp")
tc=cbind(sim.b[,c(1,3)],rep('Competitor Species',nrow(sim.b)),deparse.level = 0)
colnames(tc)=c("Time","Abund","sp")
panB=as.data.frame(rbind(tp,tc))
transTime.b=which(abs(sim.b[,2]-sim.b[,3])==min(abs(sim.b[,2]-sim.b[,3])))
transN.b=sim[transTime.b,2]

#stock sp1
tp=cbind(sim.c[,c(1,2)],rep('Preferred Species',nrow(sim.c)),deparse.level = 0)
colnames(tp)=c("Time","Abund","sp")
tc=cbind(sim.c[,c(1,3)],rep('Competitor Species',nrow(sim.c)),deparse.level = 0)
colnames(tc)=c("Time","Abund","sp")
panC=as.data.frame(rbind(tp,tc))
transTime.c=which(abs(sim.c[,2]-sim.c[,3])==min(abs(sim.c[,2]-sim.c[,3])))
transN.c=sim[transTime.c,2]

#stock and harvest
tp=cbind(sim.d[,c(1,2)],rep('Preferred Species',nrow(sim.d)),deparse.level = 0)
colnames(tp)=c("Time","Abund","sp")
tc=cbind(sim.d[,c(1,3)],rep('Competitor Species',nrow(sim.d)),deparse.level = 0)
colnames(tc)=c("Time","Abund","sp")
panD=as.data.frame(rbind(tp,tc))
transTime.d=which(abs(sim.d[,2]-sim.d[,3])==min(abs(sim.d[,2]-sim.d[,3])))
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

alt.fig4=ggarrange(a4,b4,c4,d4,labels = c("a","b","c","d"),common.legend = T,legend = "top",hjust = c(0.5,1,0.5,1))
alt.fig4=annotate_figure(alt.fig4,
                left = text_grob("Adult Abundance", rot = 90, vjust=0),
                bottom = text_grob("Time"))
alt.fig4
ggsave('alt.fig4_revision.jpeg', plot = alt.fig4, dpi = 700)
