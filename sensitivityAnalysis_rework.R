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
    comp=(temp$A2/temp$A1)<0.6 # is A2 abund <60% of A1?
    harv=c(min(temp$qEs[comp!=init]),max(temp$qEs[comp!=init])) # pulling max and min harv rates that cause a flip
    tips$parVal[i]=unique(x[,2])[i]
    tips$minQTP[i]=harv[1]
    tips$maxQTP[i]=harv[2]
    tips$dom[i]=ifelse(init==T,"A1","A2")
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
store=data.frame(A1=0,A2=0,J1=0,J2=0);store=cbind(combos,store)
store2=data.frame(A1=0,A2=0,J1=0,J2=0);store2=cbind(combos,store2)

y01=c(500,5000,0,0)
for(i in 1:nrow(store)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(0.1,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(c(s1=0.1*store$ss[i],m1=0.1,cJ1A1=0.001,cJ1A2=0.05,cJ1J2=0.003,v1=1),
      c(s2=0.1,m2=0.1,cJ2A2=0.001,cJ2A1=0.03,cJ2J1=0.003,v2=1))
  sim=ode(y=y01,times=times,func=simBiggsQ2,parms=p)
  store$A1[i]=sim[nrow(sim)-1,2]
  store$A2[i]=sim[nrow(sim)-1,3]
  store$J1[i]=sim[nrow(sim)-1,4]
  store$J2[i]=sim[nrow(sim)-1,5]
  #store$ss[i]=ss[f]
}
y02=c(5000,500,0,0)
for(i in 1:nrow(store2)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store2$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(0.1,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(c(s1=0.1*store2$ss[i],m1=0.1,cJ1A1=0.001,cJ1A2=0.05,cJ1J2=0.003,v1=1),
      c(s2=0.1,m2=0.1,cJ2A2=0.001,cJ2A1=0.03,cJ2J1=0.003,v2=1))
  sim=ode(y=y02,times=times,func=simBiggsQ2,parms=p)
  store2$A1[i]=sim[nrow(sim)-1,2]
  store2$A2[i]=sim[nrow(sim)-1,3]
  store2$J1[i]=sim[nrow(sim)-1,4]
  store2$J2[i]=sim[nrow(sim)-1,5]
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
store=data.frame(A1=0,A2=0,J1=0,J2=0);store=cbind(combos,store)
store2=data.frame(A1=0,A2=0,J1=0,J2=0);store2=cbind(combos,store2)

y01=c(500,5000,0,0)
for(i in 1:nrow(store)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(0.1,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(c(s1=0.1,m1=0.1*store$ms[i],cJ1A1=0.001,cJ1A2=0.05,cJ1J2=0.003,v1=1),
      c(s2=0.1,m2=0.1,cJ2A2=0.001,cJ2A1=0.03,cJ2J1=0.003,v2=1))
  sim=ode(y=y01,times=times,func=simBiggsQ2,parms=p)
  store$A1[i]=sim[nrow(sim)-1,2]
  store$A2[i]=sim[nrow(sim)-1,3]
  store$J1[i]=sim[nrow(sim)-1,4]
  store$J2[i]=sim[nrow(sim)-1,5]
}
y02=c(5000,500,0,0)
for(i in 1:nrow(store2)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store2$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(0.1,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(c(s1=0.1,m1=0.1*store2$ms[i],cJ1A1=0.001,cJ1A2=0.05,cJ1J2=0.003,v1=1),
      c(s2=0.1,m2=0.1,cJ2A2=0.001,cJ2A1=0.03,cJ2J1=0.003,v2=1))
  sim=ode(y=y02,times=times,func=simBiggsQ2,parms=p)
  store2$A1[i]=sim[nrow(sim)-1,2]
  store2$A2[i]=sim[nrow(sim)-1,3]
  store2$J1[i]=sim[nrow(sim)-1,4]
  store2$J2[i]=sim[nrow(sim)-1,5]
}

t3=tpS(x=store,y0=y01,parmName="ms")
t4=tpS(x=store2,y0=y02,parmName="ms")



#### ADULT PREDATION ON J2 ####
#looking to see how variation in the adult predation j2 for species 1 changes whether or not stable states occur
times=1:100

cj2a1s=seq(.5,1.5,length.out = 5)
qEs=seq(0,0.3,length.out=30)
combos=expand.grid(qEs,cj2a1s);colnames(combos)=c("qEs","cj2a1s")
store=data.frame(A1=0,A2=0,J1=0,J2=0);store=cbind(combos,store)
store2=data.frame(A1=0,A2=0,J1=0,J2=0);store2=cbind(combos,store2)

y01=c(500,5000,0,0)
for(i in 1:nrow(store)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(0.1,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(c(s1=0.1,m1=0.1,cJ1A1=0.001,cJ1A2=0.05,cJ1J2=0.003,v1=1),
      c(s2=0.1,m2=0.1,cJ2A2=0.001,cJ2A1=0.03*store$cj2a1s[i],cJ2J1=0.003,v2=1))
  sim=ode(y=y01,times=times,func=simBiggsQ2,parms=p)
  store$A1[i]=sim[nrow(sim)-1,2]
  store$A2[i]=sim[nrow(sim)-1,3]
  store$J1[i]=sim[nrow(sim)-1,4]
  store$J2[i]=sim[nrow(sim)-1,5]
}
y02=c(5000,500,0,0)
for(i in 1:nrow(store2)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store2$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(0.1,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(c(s1=0.1,m1=0.1,cJ1A1=0.001,cJ1A2=0.05,cJ1J2=0.003,v1=1),
      c(s2=0.1,m2=0.1,cJ2A2=0.001,cJ2A1=0.03*store2$cj2a1s[i],cJ2J1=0.003,v2=1))
  sim=ode(y=y02,times=times,func=simBiggsQ2,parms=p)
  store2$A1[i]=sim[nrow(sim)-1,2]
  store2$A2[i]=sim[nrow(sim)-1,3]
  store2$J1[i]=sim[nrow(sim)-1,4]
  store2$J2[i]=sim[nrow(sim)-1,5]
}

t5=tpS(x=store,y0=y01,parmName="cj2a1s")
t6=tpS(x=store2,y0=y02,parmName="cj2a1s")



#### SP1 CANNIBALISM RATE ####

times=1:100

can=seq(.5,1.5,length.out = 5)
qEs=seq(0,0.3,length.out=30)
combos=expand.grid(qEs,can);colnames(combos)=c("qEs","can")
store=data.frame(A1=0,A2=0,J1=0,J2=0);store=cbind(combos,store)
store2=data.frame(A1=0,A2=0,J1=0,J2=0);store2=cbind(combos,store2)


y01=c(500,5000,0,0)
for(i in 1:nrow(store)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(0.1,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(c(s1=0.1,m1=0.1,cJ1A1=0.001*store$can[i],cJ1A2=0.05,cJ1J2=0.003,v1=1),
      c(s2=0.1,m2=0.1,cJ2A2=0.001,cJ2A1=0.03,cJ2J1=0.003,v2=1))
  sim=ode(y=y01,times=times,func=simBiggsQ2,parms=p)
  store$A1[i]=sim[nrow(sim)-1,2]
  store$A2[i]=sim[nrow(sim)-1,3]
  store$J1[i]=sim[nrow(sim)-1,4]
  store$J2[i]=sim[nrow(sim)-1,5]
}
y02=c(5000,500,0,0)
for(i in 1:nrow(store)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store2$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(0.1,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(c(s1=0.1,m1=0.1,cJ1A1=0.001*store2$can[i],cJ1A2=0.05,cJ1J2=0.003,v1=1),
      c(s2=0.1,m2=0.1,cJ2A2=0.001,cJ2A1=0.03,cJ2J1=0.003,v2=1))
  sim=ode(y=y02,times=times,func=simBiggsQ2,parms=p)
  store2$A1[i]=sim[nrow(sim)-1,2]
  store2$A2[i]=sim[nrow(sim)-1,3]
  store2$J1[i]=sim[nrow(sim)-1,4]
  store2$J2[i]=sim[nrow(sim)-1,5]
  
}

t7=tpS(x=store,y0=y01,parmName="can")
t8=tpS(x=store2,y0=y02,parmName="can")

#### EFFECT OF J1 ON J2 ####
#looking to see how variation in the j1 effect on j2 changes whether or not stable states occur
times=1:100

cj2j1s=seq(.5,1.5,length.out = 5)
qEs=seq(0,0.3,length.out=30)
combos=expand.grid(qEs,cj2a1s);colnames(combos)=c("qEs","cj2j1s")
store=data.frame(A1=0,A2=0,J1=0,J2=0);store=cbind(combos,store)
store2=data.frame(A1=0,A2=0,J1=0,J2=0);store2=cbind(combos,store2)

y01=c(500,5000,0,0)
for(i in 1:nrow(store)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(0.1,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(c(s1=0.1,m1=0.1,cJ1A1=0.001,cJ1A2=0.05,cJ1J2=0.003,v1=1),
      c(s2=0.1,m2=0.1,cJ2A2=0.001,cJ2A1=0.03,cJ2J1=0.003*store$cj2j1s[i],v2=1))
  sim=ode(y=y01,times=times,func=simBiggsQ2,parms=p)
  store$A1[i]=sim[nrow(sim)-1,2]
  store$A2[i]=sim[nrow(sim)-1,3]
  store$J1[i]=sim[nrow(sim)-1,4]
  store$J2[i]=sim[nrow(sim)-1,5]
}
y02=c(5000,500,0,0)
for(i in 1:nrow(store2)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store2$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(0.1,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(c(s1=0.1,m1=0.1,cJ1A1=0.001,cJ1A2=0.05,cJ1J2=0.003,v1=1),
      c(s2=0.1,m2=0.1,cJ2A2=0.001,cJ2A1=0.03,cJ2J1=0.003*store2$cj2j1s[i],v2=1))
  sim=ode(y=y02,times=times,func=simBiggsQ2,parms=p)
  store2$A1[i]=sim[nrow(sim)-1,2]
  store2$A2[i]=sim[nrow(sim)-1,3]
  store2$J1[i]=sim[nrow(sim)-1,4]
  store2$J2[i]=sim[nrow(sim)-1,5]
}


t9=tpS(x=store,y0=y01,parmName="cj2j1s")
t10=tpS(x=store2,y0=y02,parmName="cj2j1s")

#### EFFECT OF J2 ON J1 ####
#looking to see how variation in the j2 effect on j1 changes whether or not stable states occur
times=1:100

cj1j2s=seq(.5,1.5,length.out = 5)
qEs=seq(0,0.3,length.out=30)
combos=expand.grid(qEs,cj2a1s);colnames(combos)=c("qEs","cj1j2s")
store=data.frame(A1=0,A2=0,J1=0,J2=0);store=cbind(combos,store)
store2=data.frame(A1=0,A2=0,J1=0,J2=0);store2=cbind(combos,store2)

y01=c(500,5000,0,0)
for(i in 1:nrow(store)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(0.1,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(c(s1=0.1,m1=0.1,cJ1A1=0.001,cJ1A2=0.05,cJ1J2=0.003*store$cj1j2s[i],v1=1),
      c(s2=0.1,m2=0.1,cJ2A2=0.001,cJ2A1=0.03,cJ2J1=0.003,v2=1))
  sim=ode(y=y01,times=times,func=simBiggsQ2,parms=p)
  store$A1[i]=sim[nrow(sim)-1,2]
  store$A2[i]=sim[nrow(sim)-1,3]
  store$J1[i]=sim[nrow(sim)-1,4]
  store$J2[i]=sim[nrow(sim)-1,5]
}
y02=c(5000,500,0,0)
for(i in 1:nrow(store2)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store2$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(0.1,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(c(s1=0.1,m1=0.1,cJ1A1=0.001,cJ1A2=0.05,cJ1J2=0.003*store$cj1j2s[i],v1=1),
      c(s2=0.1,m2=0.1,cJ2A2=0.001,cJ2A1=0.03,cJ2J1=0.003,v2=1))
  sim=ode(y=y02,times=times,func=simBiggsQ2,parms=p)
  store2$A1[i]=sim[nrow(sim)-1,2]
  store2$A2[i]=sim[nrow(sim)-1,3]
  store2$J1[i]=sim[nrow(sim)-1,4]
  store2$J2[i]=sim[nrow(sim)-1,5]
}


t11=tpS(x=store,y0=y01,parmName="cj1j2s")
t12=tpS(x=store2,y0=y02,parmName="cj1j2s")


#### ADULT PREDATION ON J1 ####
#looking to see how variation in the adult predation j1 for species 2 changes whether or not stable states occur
times=1:100

cj1a2s=seq(.5,1.5,length.out = 5)
qEs=seq(0,0.3,length.out=30)
combos=expand.grid(qEs,cj2a1s);colnames(combos)=c("qEs","cj1a2s")
store=data.frame(A1=0,A2=0,J1=0,J2=0);store=cbind(combos,store)
store2=data.frame(A1=0,A2=0,J1=0,J2=0);store2=cbind(combos,store2)

y01=c(500,5000,0,0)
for(i in 1:nrow(store)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(0.1,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(c(s1=0.1,m1=0.1,cJ1A1=0.001,cJ1A2=0.05*store$cj1a2s[i],cJ1J2=0.003,v1=1),
      c(s2=0.1,m2=0.1,cJ2A2=0.001,cJ2A1=0.03,cJ2J1=0.003,v2=1))
  sim=ode(y=y01,times=times,func=simBiggsQ2,parms=p)
  store$A1[i]=sim[nrow(sim)-1,2]
  store$A2[i]=sim[nrow(sim)-1,3]
  store$J1[i]=sim[nrow(sim)-1,4]
  store$J2[i]=sim[nrow(sim)-1,5]
}
y02=c(5000,500,0,0)
for(i in 1:nrow(store2)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store2$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(0.1,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(c(s1=0.1,m1=0.1,cJ1A1=0.001,cJ1A2=0.05*store2$cj1a2s[i],cJ1J2=0.003,v1=1),
      c(s2=0.1,m2=0.1,cJ2A2=0.001,cJ2A1=0.03,cJ2J1=0.003,v2=1))
  sim=ode(y=y02,times=times,func=simBiggsQ2,parms=p)
  store2$A1[i]=sim[nrow(sim)-1,2]
  store2$A2[i]=sim[nrow(sim)-1,3]
  store2$J1[i]=sim[nrow(sim)-1,4]
  store2$J2[i]=sim[nrow(sim)-1,5]
}


t13=tpS(x=store,y0=y01,parmName="cj1a2s")
t14=tpS(x=store2,y0=y02,parmName="cj1a2s")


#### B-H SENSITIVITY ####
#varying the  parms a & b for species 1 to see the effect of differing recruitment on stable states
# I need to create a few functions that have different parms since the function has them coded into it.
# based on Hilborn and Walters book parm a is maybe more important to vary
#species 2 stock recruitment will stay the same through all variants of speces 1 recruitment
levs=seq(.5,1.5,length.out = 5)
#original a=100; b=50
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

#var1 a=1; b=50
var1<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
    dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
    dJ1dt=-cJ1J2*J2*J1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))-((cJ1A1*v1*A1*J1)/(h1Fun(t)+v1+cJ1A1*A1))-s1*J1+(100*levs[1]*A1/(50+A1))+st1Fun(t)
    dJ2dt=-cJ2J1*J1*J2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))-((cJ2A2*v2*A2*J2)/(h2Fun(t)+v2+cJ2A2*A2))-s2*J2+(100*A2/(50+A2))+st2Fun(t)
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}

# var2 a=50.75; b=50
var2<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
    dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
    dJ1dt=-cJ1J2*J2*J1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))-((cJ1A1*v1*A1*J1)/(h1Fun(t)+v1+cJ1A1*A1))-s1*J1+(100*levs[2]*A1/(50+A1))+st1Fun(t)
    dJ2dt=-cJ2J1*J1*J2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))-((cJ2A2*v2*A2*J2)/(h2Fun(t)+v2+cJ2A2*A2))-s2*J2+(100*A2/(50+A2))+st2Fun(t)
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}

# var 3 a=100.5; b=50
var3<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
    dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
    dJ1dt=-cJ1J2*J2*J1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))-((cJ1A1*v1*A1*J1)/(h1Fun(t)+v1+cJ1A1*A1))-s1*J1+(100*levs[3]*A1/(50+A1))+st1Fun(t)
    dJ2dt=-cJ2J1*J1*J2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))-((cJ2A2*v2*A2*J2)/(h2Fun(t)+v2+cJ2A2*A2))-s2*J2+(100*A2/(50+A2))+st2Fun(t)
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}

# var 4 a=150.25; b=50
var4<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
    dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
    dJ1dt=-cJ1J2*J2*J1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))-((cJ1A1*v1*A1*J1)/(h1Fun(t)+v1+cJ1A1*A1))-s1*J1+(100*levs[4]*A1/(50+A1))+st1Fun(t)
    dJ2dt=-cJ2J1*J1*J2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))-((cJ2A2*v2*A2*J2)/(h2Fun(t)+v2+cJ2A2*A2))-s2*J2+(100*A2/(50+A2))+st2Fun(t)
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}

# var 5 a=200; b=50
var5<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
    dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
    dJ1dt=-cJ1J2*J2*J1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))-((cJ1A1*v1*A1*J1)/(h1Fun(t)+v1+cJ1A1*A1))-s1*J1+(100*levs[5]*A1/(50+A1))+st1Fun(t)
    dJ2dt=-cJ2J1*J1*J2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))-((cJ2A2*v2*A2*J2)/(h2Fun(t)+v2+cJ2A2*A2))-s2*J2+(100*A2/(50+A2))+st2Fun(t)
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}

times=1:100
funlist=list(simBiggsQ2,var1,var2,var3,var4,var5)
as=c(100,levs*100)
qEs=seq(0,0.3,length.out=30)
combos=expand.grid(qEs,as);colnames(combos)=c("qEs","as")
store=data.frame(A1=0,A2=0,J1=0,J2=0);store=cbind(combos,store)
store2=data.frame(A1=0,A2=0,J1=0,J2=0);store2=cbind(combos,store2)

y01=c(500,5000,0,0)
for(i in 1:nrow(store)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(0.1,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(c(s1=0.1,m1=0.1,cJ1A1=0.001,cJ1A2=0.05,cJ1J2=0.003,v1=1),
      c(s2=0.1,m2=0.1,cJ2A2=0.001,cJ2A1=0.03,cJ2J1=0.003,v2=1))
  if(store$as[i]==as[1]){
    sim=ode(y=y01,times=times,func=simBiggsQ2,parms=p)
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
  
  store$A1[i]=sim[nrow(sim)-1,2]
  store$A2[i]=sim[nrow(sim)-1,3]
  store$J1[i]=sim[nrow(sim)-1,4]
  store$J2[i]=sim[nrow(sim)-1,5]
}
y02=c(5000,500,0,0)
for(i in 1:nrow(store2)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store2$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(0.1,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(c(s1=0.1,m1=0.1,cJ1A1=0.001,cJ1A2=0.05,cJ1J2=0.003,v1=1),
      c(s2=0.1,m2=0.1,cJ2A2=0.001,cJ2A1=0.03,cJ2J1=0.003,v2=1))
  if(store2$as[i]==as[1]){
    sim=ode(y=y02,times=times,func=simBiggsQ2,parms=p)
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
  store2$A1[i]=sim[nrow(sim)-1,2]
  store2$A2[i]=sim[nrow(sim)-1,3]
  store2$J1[i]=sim[nrow(sim)-1,4]
  store2$J2[i]=sim[nrow(sim)-1,5]
}


t15=tpS(x=store,y0=y01,parmName="as")
t16=tpS(x=store2,y0=y02,parmName="as")


#### B-H PARM B SENSITIVITY ####
#original a=100; b=50
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

#var1 a=100; b=0.5
var1b<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
    dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
    dJ1dt=-cJ1J2*J2*J1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))-((cJ1A1*v1*A1*J1)/(h1Fun(t)+v1+cJ1A1*A1))-s1*J1+(100*A1/(50*levs[1]+A1))+st1Fun(t)
    dJ2dt=-cJ2J1*J1*J2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))-((cJ2A2*v2*A2*J2)/(h2Fun(t)+v2+cJ2A2*A2))-s2*J2+(100*A2/(50+A2))+st2Fun(t)
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}

# var2 a=100; b=25.375
var2b<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
    dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
    dJ1dt=-cJ1J2*J2*J1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))-((cJ1A1*v1*A1*J1)/(h1Fun(t)+v1+cJ1A1*A1))-s1*J1+(100*A1/(50*levs[2]+A1))+st1Fun(t)
    dJ2dt=-cJ2J1*J1*J2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))-((cJ2A2*v2*A2*J2)/(h2Fun(t)+v2+cJ2A2*A2))-s2*J2+(100*A2/(50+A2))+st2Fun(t)
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}

# var 3 a=100; b=50.25
var3b<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
    dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
    dJ1dt=-cJ1J2*J2*J1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))-((cJ1A1*v1*A1*J1)/(h1Fun(t)+v1+cJ1A1*A1))-s1*J1+(100*A1/(50*levs[3]+A1))+st1Fun(t)
    dJ2dt=-cJ2J1*J1*J2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))-((cJ2A2*v2*A2*J2)/(h2Fun(t)+v2+cJ2A2*A2))-s2*J2+(100*A2/(50+A2))+st2Fun(t)
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}

# var 4 a=100; b=75.125
var4b<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
    dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
    dJ1dt=-cJ1J2*J2*J1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))-((cJ1A1*v1*A1*J1)/(h1Fun(t)+v1+cJ1A1*A1))-s1*J1+(100*A1/(50*levs[4]+A1))+st1Fun(t)
    dJ2dt=-cJ2J1*J1*J2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))-((cJ2A2*v2*A2*J2)/(h2Fun(t)+v2+cJ2A2*A2))-s2*J2+(100*A2/(50+A2))+st2Fun(t)
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}

# var 5 a=100; b=100
var5b<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
    dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
    dJ1dt=-cJ1J2*J2*J1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))-((cJ1A1*v1*A1*J1)/(h1Fun(t)+v1+cJ1A1*A1))-s1*J1+(100*A1/(50*levs[5]+A1))+st1Fun(t)
    dJ2dt=-cJ2J1*J1*J2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))-((cJ2A2*v2*A2*J2)/(h2Fun(t)+v2+cJ2A2*A2))-s2*J2+(100*A2/(50+A2))+st2Fun(t)
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}

times=1:100
bs=c(50,50*levs)
qEs=seq(0,0.3,length.out=30)
combos=expand.grid(qEs,bs);colnames(combos)=c("qEs","bs")
store=data.frame(A1=0,A2=0,J1=0,J2=0);store=cbind(combos,store)
store2=data.frame(A1=0,A2=0,J1=0,J2=0);store2=cbind(combos,store2)

y01=c(500,5000,0,0)
for(i in 1:nrow(store)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(0.1,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(c(s1=0.1,m1=0.1,cJ1A1=0.001,cJ1A2=0.05,cJ1J2=0.003,v1=1),
      c(s2=0.1,m2=0.1,cJ2A2=0.001,cJ2A1=0.03,cJ2J1=0.003,v2=1))
  if(store$bs[i]==bs[1]){
    sim=ode(y=y01,times=times,func=simBiggsQ2,parms=p)
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
  
  store$A1[i]=sim[nrow(sim)-1,2]
  store$A2[i]=sim[nrow(sim)-1,3]
  store$J1[i]=sim[nrow(sim)-1,4]
  store$J2[i]=sim[nrow(sim)-1,5]
}
y02=c(5000,500,0,0)
for(i in 1:nrow(store2)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store2$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(0.1,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(c(s1=0.1,m1=0.1,cJ1A1=0.001,cJ1A2=0.05,cJ1J2=0.003,v1=1),
      c(s2=0.1,m2=0.1,cJ2A2=0.001,cJ2A1=0.03,cJ2J1=0.003,v2=1))
  if(store2$bs[i]==bs[1]){
    sim=ode(y=y02,times=times,func=simBiggsQ2,parms=p)
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
  store2$A1[i]=sim[nrow(sim)-1,2]
  store2$A2[i]=sim[nrow(sim)-1,3]
  store2$J1[i]=sim[nrow(sim)-1,4]
  store2$J2[i]=sim[nrow(sim)-1,5]
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

sen$dom=gsub("A1","Preferred Species",sen$dom); sen$dom=gsub("A2","Competitor Species",sen$dom)
cur=unique(sen$parName)
repl=c("sP","mP","cjCaP","cjPaP","cjCjP","cjPjC","cjPaC","a","b")
for(i in 1:9){
  sen$parName=gsub(cur[i],repl[i],sen$parName)
}

grd=ggplot(sen,aes(x=parVal*100,y=qETP))+theme_classic()+
  geom_point()+facet_grid(dom~parName)+
  geom_line(color="blue")+
  xlab("% of original Parmeter Value")+
  ylab("Harvest Rate")
grd
