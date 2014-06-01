
##################################################################
# Simulate VG as Gamma Time-Changed Brownian Motion ##############
##################################################################
VG.TimeChange.Euro=function(param,r,T,sim,S0,K,type){
  dt=T
  sigma=param[1]
  theta=param[2]
  nu=param[3]
  w=log(1-theta*nu-sigma^2*nu/2)/nu
  
  X=0
  dG=rgamma(sim,shape=dt/nu,scale=nu)
  Z=rnorm(sim)
  X=X+theta*dG+sigma*sqrt(dG)*Z
  S=S0*exp((r+w)*dt+X)
  
  payoff.disc=pmax(type*(S-K),0)*exp(-r*T)
  p.hat=mean(payoff.disc)
  p.se=sd(payoff.disc)/sqrt(sim)
  return(estimator=list(price=c(price=p.hat,se=p.se),
                        payoff=payoff.disc,
                        S=S,
                        dG=dG,
                        Z=Z,
                        X=X))
}


##################################################################
## Simulate VG as Gamma Difference of Gammas #####################
##################################################################
VG.DiffGamma.Euro=function(param,r,T,sim,S0,K,type){
  dt=T
  sigma=param[1]
  theta=param[2]
  nu=param[3]
  w=log(1-theta*nu-sigma^2*nu/2)/nu
  mu.p=VG.mu(sigma,theta,nu,1)
  mu.m=VG.mu(sigma,theta,nu,-1)
  
  X=0
  g.p=rgamma(sim,shape=T/nu,scale=nu*mu.p)
  g.m=rgamma(sim,shape=T/nu,scale=nu*mu.m)
  X=X+g.p-g.m
  S=S0*exp((r+w)*dt+X)
  
  payoff.disc=pmax(type*(S-K),0)*exp(-r*T)
  p.hat=mean(payoff.disc)
  p.se=sd(payoff.disc)/sqrt(sim)
  return(estimator=list(price=c(price=p.hat,se=p.se),
                        payoff=payoff.disc,
                        S=S,
                        X=X))
}


###################################################################
#this function calculates European option price via GBM Monte Carlo
####################################################################
GBM.Euro=function(r,sigma,T,sim,S0,K,type){
  dt=T
  eps=rnorm(sim)
  S=S0*exp((r-0.5*sigma^2)*dt+sigma*sqrt(dt)*eps)
  payoff.disc=pmax(type*(S-K),0)*exp(-r*T)
  p.hat=mean(payoff.disc)
  p.se=sd(payoff.disc)/sqrt(sim)
  return(list(price=c(price=p.hat,se=p.se),payoff=payoff.disc,S=S))
}



#Test the result
param=c(0.25,-0.15,0.5) #c(sigma,theta,nu)
S0=100
K=110
T=0.25
sim=1e5
r=0.05

VG.Euro(param,r,T,sim,S0,K,1)$price
VG.Diff.Euro(param,r,T,sim,S0,K,1)$price
GBM.Euro(r,param[1],T,sim,S0,K,1)$price
