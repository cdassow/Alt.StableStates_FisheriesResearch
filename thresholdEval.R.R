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





##FIGURE 2 #####

tstep=1:300

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
  #random % added or subtracted from the deterministic recruitment equation.
  set.seed(1)
  wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
  wC=rtriang(tstep, 0.5, 1.5, 1)
  wPFun=approxfun(x=tstep,y=wP,rule = 2)
  wCFun=approxfun(x=tstep,y=wC,rule = 2)
  
  p=c(sP=0.1,mP=0.1,cJPAP=0.002,cJPAC=0.05,cJPJC=0.003,vP=1,
      sC=0.1,mC=0.1,cJCAC=0.002,cJCAP=0.03,cJCJP=0.003,vC=1)
  y0=c(5000,500,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2w,parms=p)
  df$AP[i]=sim[nrow(sim)-1,2]
  df$AC[i]=sim[nrow(sim)-1,3]
  df$JP[i]=sim[nrow(sim)-1,4]
  df$JC[i]=sim[nrow(sim)-1,5]
}

df$diff=df$AP-df$AC
df$ratio=df$AC/df$AP
df$outcome.6=ifelse(df$ratio <= 0.6,"good", "bad") # does population fall below 60% of the other
df$outcome.5=ifelse(df$ratio <= 0.5,"good", "bad")
df$outcome.7=ifelse(df$ratio <= 0.7,"good", "bad")
df$outcome.8=ifelse(df$ratio <= 0.8,"good", "bad")

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
  #random % added or subtracted from the deterministic recruitment equation.
  set.seed(1)
  wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
  wC=rtriang(tstep, 0.5, 1.5, 1)
  wPFun=approxfun(x=tstep,y=wP,rule = 2)
  wCFun=approxfun(x=tstep,y=wC,rule = 2)
  
  p=c(sP=0.1,mP=0.1,cJPAP=0.002,cJPAC=0.05,cJPJC=0.003,vP=1,
      sC=0.1,mC=0.1,cJCAC=0.002,cJCAP=0.03,cJCJP=0.003,vC=1)
  y0=c(5000,500,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2w,parms=p)
  dfwo$AP[i]=sim[nrow(sim)-1,2]
  dfwo$AC[i]=sim[nrow(sim)-1,3]
  dfwo$JP[i]=sim[nrow(sim)-1,4]
  dfwo$JC[i]=sim[nrow(sim)-1,5]
}

dfwo$diff=dfwo$AP-dfwo$AC
dfwo$ratio=dfwo$AC/dfwo$AP
dfwo$outcome.6=ifelse(dfwo$ratio <= 0.6,"good", "bad") # does population fall below 60% of the other
dfwo$outcome.5=ifelse(dfwo$ratio <= 0.5,"good", "bad")
dfwo$outcome.7=ifelse(dfwo$ratio <= 0.7,"good", "bad")
dfwo$outcome.8=ifelse(dfwo$ratio <= 0.8,"good", "bad")

df$mod=rep("Maintain preferred, ignore competitor",nrow(df))
dfwo$mod=rep("Maintain preferred, harv competitor",nrow(dfwo))
allSen.w=rbind(df,dfwo)
allSen.w$percDiff=(abs(allSen.w$AP-allSen.w$AC)/((allSen.w$AP+allSen.w$AC)/2))*100


alt.fig2=ggplot(data=allSen.w, aes(x=allSen.w$X,y=allSen.w$Y,linetype=allSen.w$mod))+theme_classic()+
  geom_contour(aes(z=allSen.w$ratio),breaks = c(0.6), color='black', size=1.7)+
  geom_segment(aes(x=0.20, y=30, xend=0.25, yend=10), arrow = arrow(type='closed'),size=1)+
  labs(x="Preferred Species Harvest Rate", y="Preferred Species Stocking",linetype="Scenario")+
  theme(legend.position = 'none')+coord_cartesian(ylim=c(0,100), xlim=c(0,0.3))
alt.fig2

# comparison of different threshold values for determining 'dominance'
#60% presented in manuscript
alt.fig2=ggplot(data=allSen.w, aes(x=allSen.w$X,y=allSen.w$Y,linetype=allSen.w$mod))+theme_classic()+
  geom_contour(aes(z=allSen.w$ratio),breaks = c(0.6), color='black', size=1.7)+
  geom_segment(aes(x=0.20, y=30, xend=0.25, yend=10), arrow = arrow(type='closed'),size=1)+
  labs(x="Preferred Species Harvest Rate", y="Preferred Species Stocking",linetype="Scenario", title = '60% threshold')+
  theme(legend.position = 'none')+coord_cartesian(ylim=c(0,100), xlim=c(0,0.3))
#50%
alt.fig2.5=ggplot(data=allSen.w, aes(x=allSen.w$X,y=allSen.w$Y,linetype=allSen.w$mod))+theme_classic()+
  geom_contour(aes(z=allSen.w$ratio),breaks = c(0.5), color='black', size=1.7)+
  #geom_segment(aes(x=0.20, y=30, xend=0.25, yend=10), arrow = arrow(type='closed'),size=1)+
  labs(x="Preferred Species Harvest Rate", y="Preferred Species Stocking",linetype="Scenario", title = '50% threshold')+
  theme(legend.position = 'none')+coord_cartesian(ylim=c(0,100), xlim=c(0,0.3))
#70%
alt.fig2.7=ggplot(data=allSen.w, aes(x=allSen.w$X,y=allSen.w$Y,linetype=allSen.w$mod))+theme_classic()+
  geom_contour(aes(z=allSen.w$ratio),breaks = c(0.7), color='black', size=1.7)+
  #geom_segment(aes(x=0.20, y=30, xend=0.25, yend=10), arrow = arrow(type='closed'),size=1)+
  labs(x="Preferred Species Harvest Rate", y="Preferred Species Stocking",linetype="Scenario", title = '70% threshold')+
  theme(legend.position = 'none')+coord_cartesian(ylim=c(0,100), xlim=c(0,0.3))
#80%
alt.fig2.8=ggplot(data=allSen.w, aes(x=allSen.w$X,y=allSen.w$Y,linetype=allSen.w$mod))+theme_classic()+
  geom_contour(aes(z=allSen.w$ratio),breaks = c(0.8), color='black', size=1.7)+
  #geom_segment(aes(x=0.20, y=30, xend=0.25, yend=10), arrow = arrow(type='closed'),size=1)+
  labs(x="Preferred Species Harvest Rate", y="Preferred Species Stocking",linetype="Scenario", title = '80% threshold')+
  theme(legend.position = 'none')+coord_cartesian(ylim=c(0,100), xlim=c(0,0.3))
ggarrange(alt.fig2, alt.fig2.5,alt.fig2.7,alt.fig2.8)

# number of good and bad outcomes depending on what threshold you pick, doesn't change with threshold value, final abundances are always far enough apart that our conclusions are robust to threshold value choice.
table(allSen.w$outcome.6[allSen.w$mod=='Maintain preferred, ignore competitor'])
table(allSen.w$outcome.5[allSen.w$mod=='Maintain preferred, ignore competitor'])
table(allSen.w$outcome.7[allSen.w$mod=='Maintain preferred, ignore competitor'])
table(allSen.w$outcome.8[allSen.w$mod=='Maintain preferred, ignore competitor'])

table(allSen.w$outcome.6[allSen.w$mod=='Maintain preferred, harv competitor'])
table(allSen.w$outcome.5[allSen.w$mod=='Maintain preferred, harv competitor'])
table(allSen.w$outcome.7[allSen.w$mod=='Maintain preferred, harv competitor'])
table(allSen.w$outcome.8[allSen.w$mod=='Maintain preferred, harv competitor'])

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
  #random % added or subtracted from the deterministic recruitment equation.
  set.seed(1)
  wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
  wC=rtriang(tstep, 0.5, 1.5, 1)
  wPFun=approxfun(x=tstep,y=wP,rule = 2)
  wCFun=approxfun(x=tstep,y=wC,rule = 2)
  
  p=c(sP=0.1,mP=0.1,cJPAP=0.002,cJPAC=0.05,cJPJC=0.003,vP=1,
      sC=0.1,mC=0.1,cJCAC=0.002,cJCAP=0.03,cJCJP=0.003,vC=1)
  y0=c(5000,500,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2w,parms=p)
  dfT$AP[i]=sim[nrow(sim)-1,2]
  dfT$AC[i]=sim[nrow(sim)-1,3]
  dfT$JP[i]=sim[nrow(sim)-1,4]
  dfT$JC[i]=sim[nrow(sim)-1,5]
}

dfT$diff=dfT$AP-dfT$AC
dfT$ratio=dfT$AC/dfT$AP

alt.fig3=ggplot(data=dfT, aes(x=dfT$X,y=dfT$Y,linetype=as.factor(dfT$sp1H)))+theme_classic()+
  geom_contour(aes(z=dfT$ratio),breaks = c(0.6), color='black', size=1)+
  labs(x="Competitor Species Harvest Rate", y="Preferred Species Stocking",linetype="Preferred Species Harvest", title = '60%')+
  theme(legend.position = c(.75,.75), legend.key.width = unit(3, 'line'))+
  scale_linetype_manual(values = c("solid","dashed","twodash"))+
  coord_cartesian(xlim=c(0,0.15), ylim=c(0,80))
alt.fig3

alt.fig3.5=ggplot(data=dfT, aes(x=dfT$X,y=dfT$Y,linetype=as.factor(dfT$sp1H)))+theme_classic()+
  geom_contour(aes(z=dfT$ratio),breaks = c(0.5), color='black', size=1)+
  labs(x="Competitor Species Harvest Rate", y="Preferred Species Stocking",linetype="Preferred Species Harvest", title = '50%')+
  theme(legend.position = c(.75,.75), legend.key.width = unit(3, 'line'))+
  scale_linetype_manual(values = c("solid","dashed","twodash"))+
  coord_cartesian(xlim=c(0,0.15), ylim=c(0,80))

alt.fig3.7=ggplot(data=dfT, aes(x=dfT$X,y=dfT$Y,linetype=as.factor(dfT$sp1H)))+theme_classic()+
  geom_contour(aes(z=dfT$ratio),breaks = c(0.7), color='black', size=1)+
  labs(x="Competitor Species Harvest Rate", y="Preferred Species Stocking",linetype="Preferred Species Harvest", title = '70%')+
  theme(legend.position = c(.75,.75), legend.key.width = unit(3, 'line'))+
  scale_linetype_manual(values = c("solid","dashed","twodash"))+
  coord_cartesian(xlim=c(0,0.15), ylim=c(0,80))

alt.fig3.8=ggplot(data=dfT, aes(x=dfT$X,y=dfT$Y,linetype=as.factor(dfT$sp1H)))+theme_classic()+
  geom_contour(aes(z=dfT$ratio),breaks = c(0.8), color='black', size=1)+
  labs(x="Competitor Species Harvest Rate", y="Preferred Species Stocking",linetype="Preferred Species Harvest", title = '80%')+
  theme(legend.position = c(.75,.75), legend.key.width = unit(3, 'line'))+
  scale_linetype_manual(values = c("solid","dashed","twodash"))+
  coord_cartesian(xlim=c(0,0.15), ylim=c(0,80))
ggarrange(alt.fig3, alt.fig3.5, alt.fig3.7, alt.fig3.8)

dfT$outcome.6=ifelse(dfT$ratio<=0.6,'good','bad')
dfT$outcome.5=ifelse(dfT$ratio<=0.5,'good','bad')
dfT$outcome.7=ifelse(dfT$ratio<=0.7,'good','bad')
dfT$outcome.8=ifelse(dfT$ratio<=0.8,'good','bad')

table(dfT$outcome.6[dfT$sp1H==0.05])
table(dfT$outcome.5[dfT$sp1H==0.05])
table(dfT$outcome.7[dfT$sp1H==0.05])
table(dfT$outcome.8[dfT$sp1H==0.05])

table(dfT$outcome.6[dfT$sp1H==0.10])
table(dfT$outcome.5[dfT$sp1H==0.10])
table(dfT$outcome.7[dfT$sp1H==0.10])
table(dfT$outcome.8[dfT$sp1H==0.10])

table(dfT$outcome.6[dfT$sp1H==0.15])
table(dfT$outcome.5[dfT$sp1H==0.15])
table(dfT$outcome.7[dfT$sp1H==0.15])
table(dfT$outcome.8[dfT$sp1H==0.15])
