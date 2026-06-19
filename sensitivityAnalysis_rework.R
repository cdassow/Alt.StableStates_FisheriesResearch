## CJD 1.12.2021
## Updated sensitivity analysis with better plots, this is supplementary figure 1

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

#finding transition points for sensitivity analysis
tpS=function(x,y0,parmName){
  tips=data.frame(parVal=numeric(length(unique(x[,2]))),
                  minQTP=numeric(length(unique(x[,2]))),
                  maxQTP=numeric(length(unique(x[,2]))),
                  dom=character(length(unique(x[,2]))),
                  parName=rep(parmName,length(unique(x[,2]))))
  fits=data.frame(parName=parmName,intercept=numeric(1),slope=numeric(1))
  for(i in 1:length(unique(x[,2]))){
    init=y0[1]>y0[2] # initial state
    temp=x[x[,2]==unique(x[,2])[i],] # pulling rows of x with the parm value [i]
    comp=(temp$AC/temp$AP)<0.6 # is AC abund <60% of AP?
    harv=c(min(temp$qEs[comp!=init]),max(temp$qEs[comp!=init])) # pulling max and min harv rates that cause a flip
    tips$parVal[i]=unique(x[,2])[i]
    tips$minQTP[i]=harv[1]
    tips$maxQTP[i]=harv[2]
    tips$dom[i]=ifelse(init==T,"AP","AC")
  }
  if(length(unique(tips$minQTP[is.finite(tips$minQTP)]))>1){ # fit model to relationship between tipping point and parm values
    fit=lm(tips$minQTP[is.finite(tips$minQTP)]~tips$parVal[is.finite(tips$minQTP)])
  }else{fit=lm(tips$maxQTP[is.finite(tips$maxQTP)]~tips$parVal[is.finite(tips$maxQTP)])}
  fits$intercept=fit$coefficients[1]
  fits$slope=fit$coefficients[2]
  return(list(tips,fits))
}

#### JUVENILE SURVIVAL RATE ####
tstep=times=1:100
ss=seq(.5,1.5,length.out = 5)
qEs=seq(0,0.3,length.out=30)
combos=expand.grid(qEs,ss);colnames(combos)=c("qEs","ss")
store=data.frame(AP=0,AC=0,JP=0,JC=0);store=cbind(combos,store)
store2=data.frame(AP=0,AC=0,JP=0,JC=0);store2=cbind(combos,store2)

y01=c(500,5000,0,0)
for(i in 1:nrow(store)){
  qEPFun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))),rule = 2)
  qECFun=approxfun(x=tstep,y=rep(0.1,length(tstep)),rule = 2)
  v.PFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  v.CFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  kPFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  kCFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  #random % added or subtracted from the deterministic recruitment equation.
  set.seed(1)
  wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
  wC=rtriang(tstep, 0.5, 1.5, 1)
  wPFun=approxfun(x=tstep,y=wP,rule = 2)
  wCFun=approxfun(x=tstep,y=wC,rule = 2)
  
  p=c(c(sP=0.1*store$ss[i],mP=0.1,cJPAP=0.001,cJPAC=0.05,cJPJC=0.003,vP=1),
      c(sC=0.1,mC=0.1,cJCAC=0.001,cJCAP=0.03,cJCJP=0.003,vC=1))
  sim=ode(y=y01,times=times,func=simBiggsQ2w,parms=p)
  store$AP[i]=sim[nrow(sim)-1,2]
  store$AC[i]=sim[nrow(sim)-1,3]
  store$JP[i]=sim[nrow(sim)-1,4]
  store$JC[i]=sim[nrow(sim)-1,5]
  #store$ss[i]=ss[f]
}
y02=c(5000,500,0,0)
for(i in 1:nrow(store2)){
  qEPFun=approxfun(x=tstep,y=c(rep(store2$qEs[i],length(tstep))),rule = 2)
  qECFun=approxfun(x=tstep,y=rep(0.1,length(tstep)),rule = 2)
  v.PFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  v.CFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  kPFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  kCFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  #random % added or subtracted from the deterministic recruitment equation.
  set.seed(1)
  wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
  wC=rtriang(tstep, 0.5, 1.5, 1)
  wPFun=approxfun(x=tstep,y=wP,rule = 2)
  wCFun=approxfun(x=tstep,y=wC,rule = 2)
  
  p=c(c(sP=0.1*store2$ss[i],mP=0.1,cJPAP=0.001,cJPAC=0.05,cJPJC=0.003,vP=1),
      c(sC=0.1,mC=0.1,cJCAC=0.001,cJCAP=0.03,cJCJP=0.003,vC=1))
  sim=ode(y=y02,times=times,func=simBiggsQ2w,parms=p)
  store2$AP[i]=sim[nrow(sim)-1,2]
  store2$AC[i]=sim[nrow(sim)-1,3]
  store2$JP[i]=sim[nrow(sim)-1,4]
  store2$JC[i]=sim[nrow(sim)-1,5]
  #store2$ss[i]=ss[f]
}

t1=tpS(x=store,y0=y01,parmName="ss")
t2=tpS(x=store2,y0=y02,parmName="ss")




#### ADULT NATURAL MORTALITY RATE ####
#looking to see how variation in the adult natural mortality rate for species 1 changes whether or not stable states occur
times=1:100

ms=seq(.5,1.5,length.out = 5)
qEs=seq(0,0.3,length.out=30)
combos=expand.grid(qEs,ms);colnames(combos)=c("qEs","ms")
store=data.frame(AP=0,AC=0,JP=0,JC=0);store=cbind(combos,store)
store2=data.frame(AP=0,AC=0,JP=0,JC=0);store2=cbind(combos,store2)

y01=c(500,5000,0,0)
for(i in 1:nrow(store)){
  qEPFun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))),rule = 2)
  qECFun=approxfun(x=tstep,y=rep(0.1,length(tstep)),rule = 2)
  v.PFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  v.CFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  kPFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  kCFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  #random % added or subtracted from the deterministic recruitment equation.
  set.seed(1)
  wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
  wC=rtriang(tstep, 0.5, 1.5, 1)
  wPFun=approxfun(x=tstep,y=wP,rule = 2)
  wCFun=approxfun(x=tstep,y=wC,rule = 2)
  
  
  p=c(c(sP=0.1,mP=0.1*store$ms[i],cJPAP=0.001,cJPAC=0.05,cJPJC=0.003,vP=1),
      c(sC=0.1,mC=0.1,cJCAC=0.001,cJCAP=0.03,cJCJP=0.003,vC=1))
  sim=ode(y=y01,times=times,func=simBiggsQ2w,parms=p)
  store$AP[i]=sim[nrow(sim)-1,2]
  store$AC[i]=sim[nrow(sim)-1,3]
  store$JP[i]=sim[nrow(sim)-1,4]
  store$JC[i]=sim[nrow(sim)-1,5]
}
y02=c(5000,500,0,0)
for(i in 1:nrow(store2)){
  qEPFun=approxfun(x=tstep,y=c(rep(store2$qEs[i],length(tstep))),rule = 2)
  qECFun=approxfun(x=tstep,y=rep(0.1,length(tstep)),rule = 2)
  v.PFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  v.CFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  kPFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  kCFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  #random % added or subtracted from the deterministic recruitment equation.
  set.seed(1)
  wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
  wC=rtriang(tstep, 0.5, 1.5, 1)
  wPFun=approxfun(x=tstep,y=wP,rule = 2)
  wCFun=approxfun(x=tstep,y=wC,rule = 2)
  
  
  p=c(c(sP=0.1,mP=0.1*store2$ms[i],cJPAP=0.001,cJPAC=0.05,cJPJC=0.003,vP=1),
      c(sC=0.1,mC=0.1,cJCAC=0.001,cJCAP=0.03,cJCJP=0.003,vC=1))
  sim=ode(y=y02,times=times,func=simBiggsQ2w,parms=p)
  store2$AP[i]=sim[nrow(sim)-1,2]
  store2$AC[i]=sim[nrow(sim)-1,3]
  store2$JP[i]=sim[nrow(sim)-1,4]
  store2$JC[i]=sim[nrow(sim)-1,5]
}

t3=tpS(x=store,y0=y01,parmName="ms")
t4=tpS(x=store2,y0=y02,parmName="ms")



#### ADULT PREDATION ON JC ####
#looking to see how variation in the adult predation jC for preferred species changes whether or not stable states occur
times=1:100

cjCaPs=seq(.5,1.5,length.out = 5)
qEs=seq(0,0.3,length.out=30)
combos=expand.grid(qEs,cjCaPs);colnames(combos)=c("qEs","cjCaPs")
store=data.frame(AP=0,AC=0,JP=0,JC=0);store=cbind(combos,store)
store2=data.frame(AP=0,AC=0,JP=0,JC=0);store2=cbind(combos,store2)

y01=c(500,5000,0,0)
for(i in 1:nrow(store)){
  qEPFun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))),rule = 2)
  qECFun=approxfun(x=tstep,y=rep(0.1,length(tstep)),rule = 2)
  v.PFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  v.CFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  kPFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  kCFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  #random % added or subtracted from the deterministic recruitment equation.
  set.seed(1)
  wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
  wC=rtriang(tstep, 0.5, 1.5, 1)
  wPFun=approxfun(x=tstep,y=wP,rule = 2)
  wCFun=approxfun(x=tstep,y=wC,rule = 2)
  
  
  p=c(c(sP=0.1,mP=0.1,cJPAP=0.001,cJPAC=0.05,cJPJC=0.003,vP=1),
      c(sC=0.1,mC=0.1,cJCAC=0.001,cJCAP=0.03*store$cjCaPs[i],cJCJP=0.003,vC=1))
  sim=ode(y=y01,times=times,func=simBiggsQ2w,parms=p)
  store$AP[i]=sim[nrow(sim)-1,2]
  store$AC[i]=sim[nrow(sim)-1,3]
  store$JP[i]=sim[nrow(sim)-1,4]
  store$JC[i]=sim[nrow(sim)-1,5]
}
y02=c(5000,500,0,0)
for(i in 1:nrow(store2)){
  qEPFun=approxfun(x=tstep,y=c(rep(store2$qEs[i],length(tstep))),rule = 2)
  qECFun=approxfun(x=tstep,y=rep(0.1,length(tstep)),rule = 2)
  v.PFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  v.CFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  kPFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  kCFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  #random % added or subtracted from the deterministic recruitment equation.
  set.seed(1)
  wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
  wC=rtriang(tstep, 0.5, 1.5, 1)
  wPFun=approxfun(x=tstep,y=wP,rule = 2)
  wCFun=approxfun(x=tstep,y=wC,rule = 2)
  
  
  p=c(c(sP=0.1,mP=0.1,cJPAP=0.001,cJPAC=0.05,cJPJC=0.003,vP=1),
      c(sC=0.1,mC=0.1,cJCAC=0.001,cJCAP=0.03*store2$cjCaPs[i],cJCJP=0.003,vC=1))
  sim=ode(y=y02,times=times,func=simBiggsQ2w,parms=p)
  store2$AP[i]=sim[nrow(sim)-1,2]
  store2$AC[i]=sim[nrow(sim)-1,3]
  store2$JP[i]=sim[nrow(sim)-1,4]
  store2$JC[i]=sim[nrow(sim)-1,5]
}

t5=tpS(x=store,y0=y01,parmName="cjCaPs")
t6=tpS(x=store2,y0=y02,parmName="cjCaPs")



#### SP1 CANNIBALISM RATE ####

times=1:100

can=seq(.5,1.5,length.out = 5)
qEs=seq(0,0.3,length.out=30)
combos=expand.grid(qEs,can);colnames(combos)=c("qEs","can")
store=data.frame(AP=0,AC=0,JP=0,JC=0);store=cbind(combos,store)
store2=data.frame(AP=0,AC=0,JP=0,JC=0);store2=cbind(combos,store2)


y01=c(500,5000,0,0)
for(i in 1:nrow(store)){
  qEPFun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))),rule = 2)
  qECFun=approxfun(x=tstep,y=rep(0.1,length(tstep)),rule = 2)
  v.PFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  v.CFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  kPFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  kCFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  #random % added or subtracted from the deterministic recruitment equation.
  set.seed(1)
  wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
  wC=rtriang(tstep, 0.5, 1.5, 1)
  wPFun=approxfun(x=tstep,y=wP,rule = 2)
  wCFun=approxfun(x=tstep,y=wC,rule = 2)
  
  
  p=c(c(sP=0.1,mP=0.1,cJPAP=0.001*store$can[i],cJPAC=0.05,cJPJC=0.003,vP=1),
      c(sC=0.1,mC=0.1,cJCAC=0.001,cJCAP=0.03,cJCJP=0.003,vC=1))
  sim=ode(y=y01,times=times,func=simBiggsQ2w,parms=p)
  store$AP[i]=sim[nrow(sim)-1,2]
  store$AC[i]=sim[nrow(sim)-1,3]
  store$JP[i]=sim[nrow(sim)-1,4]
  store$JC[i]=sim[nrow(sim)-1,5]
}
y02=c(5000,500,0,0)
for(i in 1:nrow(store)){
  qEPFun=approxfun(x=tstep,y=c(rep(store2$qEs[i],length(tstep))),rule = 2)
  qECFun=approxfun(x=tstep,y=rep(0.1,length(tstep)),rule = 2)
  v.PFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  v.CFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  kPFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  kCFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  #random % added or subtracted from the deterministic recruitment equation.
  set.seed(1)
  wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
  wC=rtriang(tstep, 0.5, 1.5, 1)
  wPFun=approxfun(x=tstep,y=wP,rule = 2)
  wCFun=approxfun(x=tstep,y=wC,rule = 2)
  
  
  p=c(c(sP=0.1,mP=0.1,cJPAP=0.001*store2$can[i],cJPAC=0.05,cJPJC=0.003,vP=1),
      c(sC=0.1,mC=0.1,cJCAC=0.001,cJCAP=0.03,cJCJP=0.003,vC=1))
  sim=ode(y=y02,times=times,func=simBiggsQ2w,parms=p)
  store2$AP[i]=sim[nrow(sim)-1,2]
  store2$AC[i]=sim[nrow(sim)-1,3]
  store2$JP[i]=sim[nrow(sim)-1,4]
  store2$JC[i]=sim[nrow(sim)-1,5]
  
}

t7=tpS(x=store,y0=y01,parmName="can")
t8=tpS(x=store2,y0=y02,parmName="can")

#### EFFECT OF JP ON JC ####
#looking to see how variation in the jP effect on jC changes whether or not stable states occur
times=1:100

cjCjPs=seq(.5,1.5,length.out = 5)
qEs=seq(0,0.3,length.out=30)
combos=expand.grid(qEs,cjCjPs);colnames(combos)=c("qEs","cjCjPs")
store=data.frame(AP=0,AC=0,JP=0,JC=0);store=cbind(combos,store)
store2=data.frame(AP=0,AC=0,JP=0,JC=0);store2=cbind(combos,store2)

y01=c(500,5000,0,0)
for(i in 1:nrow(store)){
  qEPFun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))),rule = 2)
  qECFun=approxfun(x=tstep,y=rep(0.1,length(tstep)),rule = 2)
  v.PFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  v.CFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  kPFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  kCFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  #random % added or subtracted from the deterministic recruitment equation.
  set.seed(1)
  wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
  wC=rtriang(tstep, 0.5, 1.5, 1)
  wPFun=approxfun(x=tstep,y=wP,rule = 2)
  wCFun=approxfun(x=tstep,y=wC,rule = 2)
  
  
  p=c(c(sP=0.1,mP=0.1,cJPAP=0.001,cJPAC=0.05,cJPJC=0.003,vP=1),
      c(sC=0.1,mC=0.1,cJCAC=0.001,cJCAP=0.03,cJCJP=0.003*store$cjCjPs[i],vC=1))
  sim=ode(y=y01,times=times,func=simBiggsQ2w,parms=p)
  store$AP[i]=sim[nrow(sim)-1,2]
  store$AC[i]=sim[nrow(sim)-1,3]
  store$JP[i]=sim[nrow(sim)-1,4]
  store$JC[i]=sim[nrow(sim)-1,5]
}
y02=c(5000,500,0,0)
for(i in 1:nrow(store2)){
  qEPFun=approxfun(x=tstep,y=c(rep(store2$qEs[i],length(tstep))),rule = 2)
  qECFun=approxfun(x=tstep,y=rep(0.1,length(tstep)),rule = 2)
  v.PFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  v.CFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  kPFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  kCFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  #random % added or subtracted from the deterministic recruitment equation.
  set.seed(1)
  wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
  wC=rtriang(tstep, 0.5, 1.5, 1)
  wPFun=approxfun(x=tstep,y=wP,rule = 2)
  wCFun=approxfun(x=tstep,y=wC,rule = 2)
  
  
  p=c(c(sP=0.1,mP=0.1,cJPAP=0.001,cJPAC=0.05,cJPJC=0.003,vP=1),
      c(sC=0.1,mC=0.1,cJCAC=0.001,cJCAP=0.03,cJCJP=0.003*store2$cjCjPs[i],vC=1))
  sim=ode(y=y02,times=times,func=simBiggsQ2w,parms=p)
  store2$AP[i]=sim[nrow(sim)-1,2]
  store2$AC[i]=sim[nrow(sim)-1,3]
  store2$JP[i]=sim[nrow(sim)-1,4]
  store2$JC[i]=sim[nrow(sim)-1,5]
}


t9=tpS(x=store,y0=y01,parmName="cjCjPs")
t10=tpS(x=store2,y0=y02,parmName="cjCjPs")

#### EFFECT OF JC ON JP ####
#looking to see how variation in the jC effect on jP changes whether or not stable states occur
times=1:100

cjPjCs=seq(.5,1.5,length.out = 5)
qEs=seq(0,0.3,length.out=30)
combos=expand.grid(qEs,cjPjCs);colnames(combos)=c("qEs","cjPjCs")
store=data.frame(AP=0,AC=0,JP=0,JC=0);store=cbind(combos,store)
store2=data.frame(AP=0,AC=0,JP=0,JC=0);store2=cbind(combos,store2)

y01=c(500,5000,0,0)
for(i in 1:nrow(store)){
  qEPFun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))),rule = 2)
  qECFun=approxfun(x=tstep,y=rep(0.1,length(tstep)),rule = 2)
  v.PFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  v.CFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  kPFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  kCFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  #random % added or subtracted from the deterministic recruitment equation.
  set.seed(1)
  wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
  wC=rtriang(tstep, 0.5, 1.5, 1)
  wPFun=approxfun(x=tstep,y=wP,rule = 2)
  wCFun=approxfun(x=tstep,y=wC,rule = 2)
  
  
  p=c(c(sP=0.1,mP=0.1,cJPAP=0.001,cJPAC=0.05,cJPJC=0.003*store$cjPjCs[i],vP=1),
      c(sC=0.1,mC=0.1,cJCAC=0.001,cJCAP=0.03,cJCJP=0.003,vC=1))
  sim=ode(y=y01,times=times,func=simBiggsQ2w,parms=p)
  store$AP[i]=sim[nrow(sim)-1,2]
  store$AC[i]=sim[nrow(sim)-1,3]
  store$JP[i]=sim[nrow(sim)-1,4]
  store$JC[i]=sim[nrow(sim)-1,5]
}
y02=c(5000,500,0,0)
for(i in 1:nrow(store2)){
  qEPFun=approxfun(x=tstep,y=c(rep(store2$qEs[i],length(tstep))),rule = 2)
  qECFun=approxfun(x=tstep,y=rep(0.1,length(tstep)),rule = 2)
  v.PFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  v.CFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  kPFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  kCFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  #random % added or subtracted from the deterministic recruitment equation.
  set.seed(1)
  wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
  wC=rtriang(tstep, 0.5, 1.5, 1)
  wPFun=approxfun(x=tstep,y=wP,rule = 2)
  wCFun=approxfun(x=tstep,y=wC,rule = 2)
  
  
  p=c(c(sP=0.1,mP=0.1,cJPAP=0.001,cJPAC=0.05,cJPJC=0.003*store$cjPjCs[i],vP=1),
      c(sC=0.1,mC=0.1,cJCAC=0.001,cJCAP=0.03,cJCJP=0.003,vC=1))
  sim=ode(y=y02,times=times,func=simBiggsQ2w,parms=p)
  store2$AP[i]=sim[nrow(sim)-1,2]
  store2$AC[i]=sim[nrow(sim)-1,3]
  store2$JP[i]=sim[nrow(sim)-1,4]
  store2$JC[i]=sim[nrow(sim)-1,5]
}


t11=tpS(x=store,y0=y01,parmName="cjPjCs")
t12=tpS(x=store2,y0=y02,parmName="cjPjCs")


#### ADULT PREDATION ON JP ####
#looking to see how variation in the adult predation jP for species 2 changes whether or not stable states occur
times=1:100

cjPaCs=seq(.5,1.5,length.out = 5)
qEs=seq(0,0.3,length.out=30)
combos=expand.grid(qEs,cjPaCs);colnames(combos)=c("qEs","cjPaCs")
store=data.frame(AP=0,AC=0,JP=0,JC=0);store=cbind(combos,store)
store2=data.frame(AP=0,AC=0,JP=0,JC=0);store2=cbind(combos,store2)

y01=c(500,5000,0,0)
for(i in 1:nrow(store)){
  qEPFun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))),rule = 2)
  qECFun=approxfun(x=tstep,y=rep(0.1,length(tstep)),rule = 2)
  v.PFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  v.CFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  kPFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  kCFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  #random % added or subtracted from the deterministic recruitment equation.
  set.seed(1)
  wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
  wC=rtriang(tstep, 0.5, 1.5, 1)
  wPFun=approxfun(x=tstep,y=wP,rule = 2)
  wCFun=approxfun(x=tstep,y=wC,rule = 2)
  
  
  p=c(c(sP=0.1,mP=0.1,cJPAP=0.001,cJPAC=0.05*store$cjPaCs[i],cJPJC=0.003,vP=1),
      c(sC=0.1,mC=0.1,cJCAC=0.001,cJCAP=0.03,cJCJP=0.003,vC=1))
  sim=ode(y=y01,times=times,func=simBiggsQ2w,parms=p)
  store$AP[i]=sim[nrow(sim)-1,2]
  store$AC[i]=sim[nrow(sim)-1,3]
  store$JP[i]=sim[nrow(sim)-1,4]
  store$JC[i]=sim[nrow(sim)-1,5]
}
y02=c(5000,500,0,0)
for(i in 1:nrow(store2)){
  qEPFun=approxfun(x=tstep,y=c(rep(store2$qEs[i],length(tstep))),rule = 2)
  qECFun=approxfun(x=tstep,y=rep(0.1,length(tstep)),rule = 2)
  v.PFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  v.CFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  kPFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  kCFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  #random % added or subtracted from the deterministic recruitment equation.
  set.seed(1)
  wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
  wC=rtriang(tstep, 0.5, 1.5, 1)
  wPFun=approxfun(x=tstep,y=wP,rule = 2)
  wCFun=approxfun(x=tstep,y=wC,rule = 2)
  
  
  p=c(c(sP=0.1,mP=0.1,cJPAP=0.001,cJPAC=0.05*store2$cjPaCs[i],cJPJC=0.003,vP=1),
      c(sC=0.1,mC=0.1,cJCAC=0.001,cJCAP=0.03,cJCJP=0.003,vC=1))
  sim=ode(y=y02,times=times,func=simBiggsQ2w,parms=p)
  store2$AP[i]=sim[nrow(sim)-1,2]
  store2$AC[i]=sim[nrow(sim)-1,3]
  store2$JP[i]=sim[nrow(sim)-1,4]
  store2$JC[i]=sim[nrow(sim)-1,5]
}


t13=tpS(x=store,y0=y01,parmName="cjPaCs")
t14=tpS(x=store2,y0=y02,parmName="cjPaCs")


#### B-H SENSITIVITY ####
#varying the  parms a & b for species 1 to see the effect of differing recruitment on stable states
# I need to create a few functions that have different parms since the function has them coded into it.
# based on Hilborn and Walters book parm a is maybe more important to vary
#species 2 stock recruitment will stay the same through all variants of speces 1 recruitment
levs=seq(.5,1.5,length.out = 5)
#original a=100; b=50
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

#var1 a=1; b=50
var1<-function(t,y,params){
  AP<-y[1]
  AC<-y[2]
  JP<-y[3]
  JC<-y[4]
  with(as.list(params),{
    dAPdt=-qEPFun(t)*AP-mP*AP+sP*JP
    dACdt=-qECFun(t)*AC-mC*AC+sC*JC
    dJPdt=-cJPJC*JC*JP-((cJPAC*vP*AC*JP)/(v.PFun(t)+vP+cJPAC*AC))-((cJPAP*vP*AP*JP)/(v.PFun(t)+vP+cJPAP*AP))-sP*JP+((100*levs[1]*AP/(50+AP))*wPFun(t))+kPFun(t)
    dJCdt=-cJCJP*JP*JC-((cJCAP*vC*AP*JC)/(v.CFun(t)+vC+cJCAP*AP))-((cJCAC*vC*AC*JC)/(v.CFun(t)+vC+cJCAC*AC))-sC*JC+((100*AC/(50+AC))*wCFun(t))+kCFun(t)
    return(list(c(dAPdt,dACdt,dJPdt,dJCdt)))
  })
}

# var2 a=50.75; b=50
var2<-function(t,y,params){
  AP<-y[1]
  AC<-y[2]
  JP<-y[3]
  JC<-y[4]
  with(as.list(params),{
    dAPdt=-qEPFun(t)*AP-mP*AP+sP*JP
    dACdt=-qECFun(t)*AC-mC*AC+sC*JC
    dJPdt=-cJPJC*JC*JP-((cJPAC*vP*AC*JP)/(v.PFun(t)+vP+cJPAC*AC))-((cJPAP*vP*AP*JP)/(v.PFun(t)+vP+cJPAP*AP))-sP*JP+((100*levs[2]*AP/(50+AP))*wPFun(t))+kPFun(t)
    dJCdt=-cJCJP*JP*JC-((cJCAP*vC*AP*JC)/(v.CFun(t)+vC+cJCAP*AP))-((cJCAC*vC*AC*JC)/(v.CFun(t)+vC+cJCAC*AC))-sC*JC+((100*AC/(50+AC))*wCFun(t))+kCFun(t)
    return(list(c(dAPdt,dACdt,dJPdt,dJCdt)))
  })
}

# var 3 a=100.5; b=50
var3<-function(t,y,params){
  AP<-y[1]
  AC<-y[2]
  JP<-y[3]
  JC<-y[4]
  with(as.list(params),{
    dAPdt=-qEPFun(t)*AP-mP*AP+sP*JP
    dACdt=-qECFun(t)*AC-mC*AC+sC*JC
    dJPdt=-cJPJC*JC*JP-((cJPAC*vP*AC*JP)/(v.PFun(t)+vP+cJPAC*AC))-((cJPAP*vP*AP*JP)/(v.PFun(t)+vP+cJPAP*AP))-sP*JP+((100*levs[3]*AP/(50+AP))*wPFun(t))+kPFun(t)
    dJCdt=-cJCJP*JP*JC-((cJCAP*vC*AP*JC)/(v.CFun(t)+vC+cJCAP*AP))-((cJCAC*vC*AC*JC)/(v.CFun(t)+vC+cJCAC*AC))-sC*JC+((100*AC/(50+AC))*wCFun(t))+kCFun(t)
    return(list(c(dAPdt,dACdt,dJPdt,dJCdt)))
  })
}

# var 4 a=150.25; b=50
var4<-function(t,y,params){
  AP<-y[1]
  AC<-y[2]
  JP<-y[3]
  JC<-y[4]
  with(as.list(params),{
    dAPdt=-qEPFun(t)*AP-mP*AP+sP*JP
    dACdt=-qECFun(t)*AC-mC*AC+sC*JC
    dJPdt=-cJPJC*JC*JP-((cJPAC*vP*AC*JP)/(v.PFun(t)+vP+cJPAC*AC))-((cJPAP*vP*AP*JP)/(v.PFun(t)+vP+cJPAP*AP))-sP*JP+((100*levs[4]*AP/(50+AP))*wPFun(t))+kPFun(t)
    dJCdt=-cJCJP*JP*JC-((cJCAP*vC*AP*JC)/(v.CFun(t)+vC+cJCAP*AP))-((cJCAC*vC*AC*JC)/(v.CFun(t)+vC+cJCAC*AC))-sC*JC+((100*AC/(50+AC))*wCFun(t))+kCFun(t)
    return(list(c(dAPdt,dACdt,dJPdt,dJCdt)))
  })
}

# var 5 a=200; b=50
var5<-function(t,y,params){
  AP<-y[1]
  AC<-y[2]
  JP<-y[3]
  JC<-y[4]
  with(as.list(params),{
    dAPdt=-qEPFun(t)*AP-mP*AP+sP*JP
    dACdt=-qECFun(t)*AC-mC*AC+sC*JC
    dJPdt=-cJPJC*JC*JP-((cJPAC*vP*AC*JP)/(v.PFun(t)+vP+cJPAC*AC))-((cJPAP*vP*AP*JP)/(v.PFun(t)+vP+cJPAP*AP))-sP*JP+((100*levs[5]*AP/(50+AP))*wPFun(t))+kPFun(t)
    dJCdt=-cJCJP*JP*JC-((cJCAP*vC*AP*JC)/(v.CFun(t)+vC+cJCAP*AP))-((cJCAC*vC*AC*JC)/(v.CFun(t)+vC+cJCAC*AC))-sC*JC+((100*AC/(50+AC))*wCFun(t))+kCFun(t)
    return(list(c(dAPdt,dACdt,dJPdt,dJCdt)))
  })
}

times=1:100
funlist=list(simBiggsQ2w,var1,var2,var3,var4,var5)
as=c(100,levs*100)
qEs=seq(0,0.3,length.out=30)
combos=expand.grid(qEs,as);colnames(combos)=c("qEs","as")
store=data.frame(AP=0,AC=0,JP=0,JC=0);store=cbind(combos,store)
store2=data.frame(AP=0,AC=0,JP=0,JC=0);store2=cbind(combos,store2)

y01=c(500,5000,0,0)
for(i in 1:nrow(store)){
  qEPFun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))),rule = 2)
  qECFun=approxfun(x=tstep,y=rep(0.1,length(tstep)),rule = 2)
  v.PFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  v.CFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  kPFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  kCFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  #random % added or subtracted from the deterministic recruitment equation.
  set.seed(1)
  wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
  wC=rtriang(tstep, 0.5, 1.5, 1)
  wPFun=approxfun(x=tstep,y=wP,rule = 2)
  wCFun=approxfun(x=tstep,y=wC,rule = 2)
  
  p=c(c(sP=0.1,mP=0.1,cJPAP=0.001,cJPAC=0.05,cJPJC=0.003,vP=1),
      c(sC=0.1,mC=0.1,cJCAC=0.001,cJCAP=0.03,cJCJP=0.003,vC=1))
  if(store$as[i]==as[1]){
    sim=ode(y=y01,times=times,func=simBiggsQ2w,parms=p)
  }else{if(store$as[i]==as[2]){
    sim=ode(y=y01,times=times,func=var5,parms=p)
  }else{if(store$as[i]==as[3]){
    sim=ode(y=y01,times=times,func=var4,parms=p)
  }else{if(store$as[i]==as[4]){
    sim=ode(y=y01,times=times,func=var3,parms=p)
  }else{if(store$as[i]==as[5]){
    sim=ode(y=y01,times=times,func=var2,parms=p)
  }else{if(store$as[i]==as[6]){
    sim=ode(y=y01,times=times,func=var1,parms=p)
  }
  }}}}}
  
  store$AP[i]=sim[nrow(sim)-1,2]
  store$AC[i]=sim[nrow(sim)-1,3]
  store$JP[i]=sim[nrow(sim)-1,4]
  store$JC[i]=sim[nrow(sim)-1,5]
}
y02=c(5000,500,0,0)
for(i in 1:nrow(store2)){
  qEPFun=approxfun(x=tstep,y=c(rep(store2$qEs[i],length(tstep))),rule = 2)
  qECFun=approxfun(x=tstep,y=rep(0.1,length(tstep)),rule = 2)
  v.PFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  v.CFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  kPFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  kCFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  #random % added or subtracted from the deterministic recruitment equation.
  set.seed(1)
  wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
  wC=rtriang(tstep, 0.5, 1.5, 1)
  wPFun=approxfun(x=tstep,y=wP,rule = 2)
  wCFun=approxfun(x=tstep,y=wC,rule = 2)
  
  p=c(c(sP=0.1,mP=0.1,cJPAP=0.001,cJPAC=0.05,cJPJC=0.003,vP=1),
      c(sC=0.1,mC=0.1,cJCAC=0.001,cJCAP=0.03,cJCJP=0.003,vC=1))
  if(store2$as[i]==as[1]){
    sim=ode(y=y02,times=times,func=simBiggsQ2w,parms=p)
  }else{if(store2$as[i]==as[2]){
    sim=ode(y=y02,times=times,func=var5,parms=p)
  }else{if(store2$as[i]==as[3]){
    sim=ode(y=y02,times=times,func=var4,parms=p)
  }else{if(store2$as[i]==as[4]){
    sim=ode(y=y02,times=times,func=var3,parms=p)
  }else{if(store2$as[i]==as[5]){
    sim=ode(y=y02,times=times,func=var2,parms=p)
  }else{if(store2$as[i]==as[6]){
    sim=ode(y=y02,times=times,func=var1,parms=p)
  }
  }}}}}
  store2$AP[i]=sim[nrow(sim)-1,2]
  store2$AC[i]=sim[nrow(sim)-1,3]
  store2$JP[i]=sim[nrow(sim)-1,4]
  store2$JC[i]=sim[nrow(sim)-1,5]
}


t15=tpS(x=store,y0=y01,parmName="as")
t16=tpS(x=store2,y0=y02,parmName="as")


#### B-H PARM B SENSITIVITY ####
#original a=100; b=50
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

#var1 a=100; b=0.5
var1b<-function(t,y,params){
  AP<-y[1]
  AC<-y[2]
  JP<-y[3]
  JC<-y[4]
  with(as.list(params),{
    dAPdt=-qEPFun(t)*AP-mP*AP+sP*JP
    dACdt=-qECFun(t)*AC-mC*AC+sC*JC
    dJPdt=-cJPJC*JC*JP-((cJPAC*vP*AC*JP)/(v.PFun(t)+vP+cJPAC*AC))-((cJPAP*vP*AP*JP)/(v.PFun(t)+vP+cJPAP*AP))-sP*JP+((100*AP/(50*levs[1]+AP))*wPFun(t))+kPFun(t)
    dJCdt=-cJCJP*JP*JC-((cJCAP*vC*AP*JC)/(v.CFun(t)+vC+cJCAP*AP))-((cJCAC*vC*AC*JC)/(v.CFun(t)+vC+cJCAC*AC))-sC*JC+((100*AC/(50+AC))*wCFun(t))+kCFun(t)
    return(list(c(dAPdt,dACdt,dJPdt,dJCdt)))
  })
}

# var2 a=100; b=25.375
var2b<-function(t,y,params){
  AP<-y[1]
  AC<-y[2]
  JP<-y[3]
  JC<-y[4]
  with(as.list(params),{
    dAPdt=-qEPFun(t)*AP-mP*AP+sP*JP
    dACdt=-qECFun(t)*AC-mC*AC+sC*JC
    dJPdt=-cJPJC*JC*JP-((cJPAC*vP*AC*JP)/(v.PFun(t)+vP+cJPAC*AC))-((cJPAP*vP*AP*JP)/(v.PFun(t)+vP+cJPAP*AP))-sP*JP+((100*AP/(50*levs[2]+AP))*wPFun(t))+kPFun(t)
    dJCdt=-cJCJP*JP*JC-((cJCAP*vC*AP*JC)/(v.CFun(t)+vC+cJCAP*AP))-((cJCAC*vC*AC*JC)/(v.CFun(t)+vC+cJCAC*AC))-sC*JC+((100*AC/(50+AC))*wCFun(t))+kCFun(t)
    return(list(c(dAPdt,dACdt,dJPdt,dJCdt)))
  })
}

# var 3 a=100; b=50.25
var3b<-function(t,y,params){
  AP<-y[1]
  AC<-y[2]
  JP<-y[3]
  JC<-y[4]
  with(as.list(params),{
    dAPdt=-qEPFun(t)*AP-mP*AP+sP*JP
    dACdt=-qECFun(t)*AC-mC*AC+sC*JC
    dJPdt=-cJPJC*JC*JP-((cJPAC*vP*AC*JP)/(v.PFun(t)+vP+cJPAC*AC))-((cJPAP*vP*AP*JP)/(v.PFun(t)+vP+cJPAP*AP))-sP*JP+((100*AP/(50*levs[3]+AP))*wPFun(t))+kPFun(t)
    dJCdt=-cJCJP*JP*JC-((cJCAP*vC*AP*JC)/(v.CFun(t)+vC+cJCAP*AP))-((cJCAC*vC*AC*JC)/(v.CFun(t)+vC+cJCAC*AC))-sC*JC+((100*AC/(50+AC))*wCFun(t))+kCFun(t)
    return(list(c(dAPdt,dACdt,dJPdt,dJCdt)))
  })
}

# var 4 a=100; b=75.125
var4b<-function(t,y,params){
  AP<-y[1]
  AC<-y[2]
  JP<-y[3]
  JC<-y[4]
  with(as.list(params),{
    dAPdt=-qEPFun(t)*AP-mP*AP+sP*JP
    dACdt=-qECFun(t)*AC-mC*AC+sC*JC
    dJPdt=-cJPJC*JC*JP-((cJPAC*vP*AC*JP)/(v.PFun(t)+vP+cJPAC*AC))-((cJPAP*vP*AP*JP)/(v.PFun(t)+vP+cJPAP*AP))-sP*JP+((100*AP/(50*levs[4]+AP))*wPFun(t))+kPFun(t)
    dJCdt=-cJCJP*JP*JC-((cJCAP*vC*AP*JC)/(v.CFun(t)+vC+cJCAP*AP))-((cJCAC*vC*AC*JC)/(v.CFun(t)+vC+cJCAC*AC))-sC*JC+((100*AC/(50+AC))*wCFun(t))+kCFun(t)
    return(list(c(dAPdt,dACdt,dJPdt,dJCdt)))
  })
}

# var 5 a=100; b=100
var5b<-function(t,y,params){
  AP<-y[1]
  AC<-y[2]
  JP<-y[3]
  JC<-y[4]
  with(as.list(params),{
    dAPdt=-qEPFun(t)*AP-mP*AP+sP*JP
    dACdt=-qECFun(t)*AC-mC*AC+sC*JC
    dJPdt=-cJPJC*JC*JP-((cJPAC*vP*AC*JP)/(v.PFun(t)+vP+cJPAC*AC))-((cJPAP*vP*AP*JP)/(v.PFun(t)+vP+cJPAP*AP))-sP*JP+((100*AP/(50*levs[5]+AP))*wPFun(t))+kPFun(t)
    dJCdt=-cJCJP*JP*JC-((cJCAP*vC*AP*JC)/(v.CFun(t)+vC+cJCAP*AP))-((cJCAC*vC*AC*JC)/(v.CFun(t)+vC+cJCAC*AC))-sC*JC+((100*AC/(50+AC))*wCFun(t))+kCFun(t)
    return(list(c(dAPdt,dACdt,dJPdt,dJCdt)))
  })
}

times=1:100
bs=c(50,50*levs)
qEs=seq(0,0.3,length.out=30)
combos=expand.grid(qEs,bs);colnames(combos)=c("qEs","bs")
store=data.frame(AP=0,AC=0,JP=0,JC=0);store=cbind(combos,store)
store2=data.frame(AP=0,AC=0,JP=0,JC=0);store2=cbind(combos,store2)

y01=c(500,5000,0,0)
for(i in 1:nrow(store)){
  qEPFun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))),rule = 2)
  qECFun=approxfun(x=tstep,y=rep(0.1,length(tstep)),rule = 2)
  v.PFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  v.CFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  kPFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  kCFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  #random % added or subtracted from the deterministic recruitment equation.
  set.seed(1)
  wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
  wC=rtriang(tstep, 0.5, 1.5, 1)
  wPFun=approxfun(x=tstep,y=wP,rule = 2)
  wCFun=approxfun(x=tstep,y=wC,rule = 2)
  
  p=c(c(sP=0.1,mP=0.1,cJPAP=0.001,cJPAC=0.05,cJPJC=0.003,vP=1),
      c(sC=0.1,mC=0.1,cJCAC=0.001,cJCAP=0.03,cJCJP=0.003,vC=1))
  if(store$bs[i]==bs[1]){
    sim=ode(y=y01,times=times,func=simBiggsQ2w,parms=p)
  }else{if(store$bs[i]==bs[2]){
    sim=ode(y=y01,times=times,func=var5b,parms=p)
  }else{if(store$bs[i]==bs[3]){
    sim=ode(y=y01,times=times,func=var4b,parms=p)
  }else{if(store$bs[i]==bs[4]){
    sim=ode(y=y01,times=times,func=var3b,parms=p)
  }else{if(store$bs[i]==bs[5]){
    sim=ode(y=y01,times=times,func=var2b,parms=p)
  }else{if(store$bs[i]==bs[6]){
    sim=ode(y=y01,times=times,func=var1b,parms=p)
  }
  }}}}}
  
  store$AP[i]=sim[nrow(sim)-1,2]
  store$AC[i]=sim[nrow(sim)-1,3]
  store$JP[i]=sim[nrow(sim)-1,4]
  store$JC[i]=sim[nrow(sim)-1,5]
}
y02=c(5000,500,0,0)
for(i in 1:nrow(store2)){
  qEPFun=approxfun(x=tstep,y=c(rep(store2$qEs[i],length(tstep))),rule = 2)
  qECFun=approxfun(x=tstep,y=rep(0.1,length(tstep)),rule = 2)
  v.PFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  v.CFun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
  kPFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  kCFun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  #random % added or subtracted from the deterministic recruitment equation.
  set.seed(1)
  wP=rtriang(tstep, 0.5, 1.5, 1) # probability distribution between 50% reduction/addition, or no change
  wC=rtriang(tstep, 0.5, 1.5, 1)
  wPFun=approxfun(x=tstep,y=wP,rule = 2)
  wCFun=approxfun(x=tstep,y=wC,rule = 2)
  
  p=c(c(sP=0.1,mP=0.1,cJPAP=0.001,cJPAC=0.05,cJPJC=0.003,vP=1),
      c(sC=0.1,mC=0.1,cJCAC=0.001,cJCAP=0.03,cJCJP=0.003,vC=1))
  if(store2$bs[i]==bs[1]){
    sim=ode(y=y02,times=times,func=simBiggsQ2w,parms=p)
  }else{if(store2$bs[i]==bs[2]){
    sim=ode(y=y02,times=times,func=var5b,parms=p)
  }else{if(store2$bs[i]==bs[3]){
    sim=ode(y=y02,times=times,func=var4b,parms=p)
  }else{if(store2$bs[i]==bs[4]){
    sim=ode(y=y02,times=times,func=var3b,parms=p)
  }else{if(store2$bs[i]==bs[5]){
    sim=ode(y=y02,times=times,func=var2b,parms=p)
  }else{if(store2$bs[i]==bs[6]){
    sim=ode(y=y02,times=times,func=var1b,parms=p)
  }
  }}}}}
  store2$AP[i]=sim[nrow(sim)-1,2]
  store2$AC[i]=sim[nrow(sim)-1,3]
  store2$JP[i]=sim[nrow(sim)-1,4]
  store2$JC[i]=sim[nrow(sim)-1,5]
}


t17=tpS(x=store,y0=y01,parmName="bs")
t18=tpS(x=store2,y0=y02,parmName="bs")

#### PLOTTING ALL RESULTS ####

sen=rbind(t1[[1]],t2[[1]],t3[[1]],t4[[1]],t5[[1]],t6[[1]],t7[[1]],t8[[1]],t9[[1]],t10[[1]],t11[[1]],t12[[1]],t13[[1]],t14[[1]],t15[[1]],t16[[1]],t17[[1]],t18[[1]])

for(i in 1:nrow(sen)){
  sen$qETP[i]=ifelse(sen$minQTP[i]==0,sen$maxQTP[i],sen$minQTP[i])
}
for(i in 1:nrow(sen)){if(sen$parVal[i]%in%as){
  if(sen$parVal[i]==as[2]){sen$parVal[i]=levs[1]}else{
    if(sen$parVal[i]==as[1]){sen$parVal[i]=1}else{
      if(sen$parVal[i]==as[3]){sen$parVal[i]=levs[2]}else{
        if(sen$parVal[i]==as[4]){sen$parVal[i]=levs[3]}else{
          if(sen$parVal[i]==as[5]){sen$parVal[i]=levs[4]}else{
            if(sen$parVal[i]==as[6]){sen$parVal[i]=levs[5]}
          }
        }
      }
    }
  }
}else{if(sen$parVal[i]==bs[2]){sen$parVal[i]=levs[1]}else{
  if(sen$parVal[i]==bs[1]){sen$parVal[i]=1}else{
    if(sen$parVal[i]==bs[3]){sen$parVal[i]=levs[2]}else{
      if(sen$parVal[i]==bs[4]){sen$parVal[i]=levs[3]}else{
        if(sen$parVal[i]==bs[5]){sen$parVal[i]=levs[4]}else{
          if(sen$parVal[i]==bs[6]){sen$parVal[i]=levs[5]}
        }
      }
    }
  }
}}}

fits=rbind(t1[[2]],t2[[2]],t3[[2]],t4[[2]],t5[[2]],t6[[2]],t7[[2]],t8[[2]],t9[[2]],t10[[2]],t11[[2]],t12[[2]],t13[[2]],t14[[2]],t15[[2]],t16[[2]],t17[[2]],t18[[2]])
dumDat=data.frame()

lms=ggplot(data = dumDat)+theme_classic()+
  geom_point()+xlim(0.01,2)+ylim(0,0.5)+
  geom_abline(data = fits, aes(slope=slope,intercept=intercept,color=parName))
lms

sen$dom=gsub("AP","Preferred Species",sen$dom); sen$dom=gsub("AC","Competitor Species",sen$dom)
cur=unique(sen$parName)
repl=c("sP","mP","cjCaP","cjPaP","cjCjP","cjPjC","cjPaC","a","b")
for(i in 1:9){
  sen$parName=gsub(cur[i],repl[i],sen$parName)
}

# rows of sen with Inf or -Inf for minQTP and maxQTP columns indicate situations where the system never flipped from competitor species to preferred species. On the plot below this results in the spike in harvest rate that goes to the top of the plot for parameters a, cJCJP, mP, and sP in the top row.

# this is noted in the supplemental information

si.fig1=ggplot(sen,aes(x=parVal*100,y=qETP))+theme_classic()+
  geom_point()+facet_grid(dom~parName)+
  geom_line(color="blue")+
  xlab("% of original Parmeter Value")+
  ylab("Harvest Rate")
si.fig1
ggsave('supplementFig1_revision.png', plot = si.fig1, dpi='print', width = 11, height = 5, units = 'in')
