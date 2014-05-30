library(VarianceGamma)
library(quantmod)

#download underlying asset prices and calculate its log returns
symbol="SPY"
dat=getSymbols(Symbols=symbol,src="yahoo")
data=SPY["2012/2014"][,6]
SPreturn=log(as.numeric(data[-1])/as.numeric(data[-length(data)]))

#estimate parameters for normal distribution
mu.norm=mean(SPreturn)
sigma.norm=sd(SPreturn)
#estimate robust parameters for normal distribution
mu.norm.rob=median(SPreturn)
sigma.norm.rob=mad(SPreturn)

#estimate parameters for VG distribution
#fit=vgFit(SPreturn,startValues="MoM")
fit=vgFit(SPreturn,startValues="US",paramStart=c(0.0009,0.006,0,0.5))
param=fit$param

#compare the density curves of different distributions
png(file="density plots.png",width=7,height=5,units="in",res=300)
plot(density(SPreturn,adjust=0.8),main="Empirical density vs. Three distributions")
curve(dnorm(x,mu.norm,sigma.norm),add=TRUE,col="blue",lty=2)
curve(dnorm(x,mu.norm.rob,sigma.norm.rob),add=TRUE,col="green",lty=3)
curve(dvg(x,param=param),add=TRUE,col="red")
legend("topright",legend=c("Empirical density","Normal density",
                           "Robust normal density","VG density"),
       col=c("black","blue","green","red"),lty=c(1,2,3,1),bty="n")
dev.off()


#compare the qqplots for normal distribution and VG distribution
#plot(fit,which=3)
png(file="QQ plots.png",width=10,height=5,units="in",res=300)
par(mfrow=c(1,2))
qqnorm(SPreturn)
qqline(SPreturn)
qqvg(SPreturn,param=param)
par(mfrow=c(1,1))
dev.off()


#this function calculates European option price via VG Monte Carlo
VG.sim=function(param,r,T,sim,S0,K,type){
  sigma=param[2]
  theta=param[3]
  nu=param[4]
  w=log(1-theta*nu-sigma^2*nu/2)/nu
  eps=rvg(sim,param=param)
  S=S0*exp((r+w)*T+eps)
  
  payoff.disc=pmax(type*(S-K),0)*exp(-r*T)
  p.hat=mean(payoff.disc)
  p.se=sd(payoff.disc)/sqrt(sim)
  return(c(price=p.hat,se=p.se))
}


#this function calculates European option price via GM Monte Carlo
GM.sim=function(r,sigma,T,sim,S0,K,type){
  eps=rnorm(sim)
  S=S0*exp((r-0.5*sigma^2)*T+sigma*sqrt(T)*eps)
  payoff.disc=pmax(type*(S-K),0)*exp(-r*T)
  p.hat=mean(payoff.disc)
  p.se=sd(payoff.disc)/sqrt(sim)
  return(c(price=p.hat,se=p.se))
}


S0=192.37
K=95
T=0.56
r=0.04
sim=100000
annual.sigma=sqrt(252)*sigma.norm

#compute the European call option price with VG and GM model
VG.sim(param,r,T,sim,S0,K,1)
GM.sim(r,annual.sigma,T,sim,S0,K,1)

