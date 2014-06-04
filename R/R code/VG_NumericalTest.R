# Initiate parameters
param=c(0.25,-0.05,0.5) #c(sigma,theta,nu)
S0=100
K=110
T=1
sim=1e5
r=0.05
steps=252
tau <- T/steps


## 1. Asian/European price for VG_MC, VG_FFT, GBM, BS
########################################################
VG.TimeChange.Euro(param,r,T,sim,S0,K,1)$price
VG.DiffGamma.Euro(param,r,T,sim,S0,K,1)$price
FFT.price(char.VG, S0 = S0, K = K, r = r, T = T, type = 1)
GBM.Euro(r,param[1],T,sim,S0,K,1)$price
bs.price(S0, K, r, T, param[1], 1) 

VG.TimeChange.Euro(param,r,T,sim,S0,K,-1)$price
VG.DiffGamma.Euro(param,r,T,sim,S0,K,-1)$price
FFT.price(char.VG, S0 = S0, K = K, r = r, T = T, type = -1)
GBM.Euro(r,param[1],T,sim,S0,K,-1)$price
bs.price(S0, K, r, T, param[1], -1)

vg.geom.asian.mc.price(S0, K, r, T, param[1], type = 1, param[3], param[2], sim, steps)$price
FFT.price.geom.asian(char.VG, S0 = S0, K = K, r = r, T = T,type = 1, steps)
bs.geom.asian.mc.price(S0, K, r, T, param[1], type = 1, sim, steps)$price
FFT.price.geom.asian(char.BS, S0 = S0, K = K, r = r, T = T,type = 1, steps)

vg.geom.asian.mc.price(S0, K, r, T, param[1], type = -1, param[3], param[2], sim, steps)$price
FFT.price.geom.asian(char.VG, S0 = S0, K = K, r = r, T = T,type = -1, steps)
bs.geom.asian.mc.price(S0, K, r, T, param[1], type = -1, sim, steps)$price
FFT.price.geom.asian(char.BS, S0 = S0, K = K, r = r, T = T,type = -1, steps)


## 2 price paths for VG and GBM
#################################################################################################
S.VG=vg.geom.asian.mc.price(S0, K, r, T, param[1], type = 1, param[3], param[2], 10, steps)$S
S.GBM=bs.geom.asian.mc.price(S0, K, r, T, param[1], type = -1, 10, steps)$S

png(file="Stock price paths for VG and GBM process.png",width=7,height=9,units="in",res=300)
par(mfrow=c(2,1))
matplot(t(S.VG),type="l",xlab="Time",ylab="Stock Price",
        main="Stock price paths for time-changed VG process
        (sigma=0.25, theta=-0.05, nu=0.5)")
matplot(t(S.GBM),type="l",xlab="Time",ylab="Stock Price",main="Stock price paths for GBM process")
dev.off()


#3 European call/put option prices change with K
###############################################################################################3
strike=seq(K-20,K+20,by=5)
p.VG=sapply(strike, function(x) VG.TimeChange.Euro(param,r,T,sim,S0,x,1)$price[1]) 
p.GBM=sapply(strike, function(x) GBM.Euro(r,param[1],T,sim,S0,x,1)$price[1])
png(file="Strike price vs call option price.png",width=7,height=5,units="in",res=300)
plot(strike,p.VG,type="l",xlab="K",ylab="Price",main="Strike price vs. Call option price")
lines(strike,p.GBM,lty=2,col="red")
legend("topright",legend=c("VG","GBM"),lty=c(1,2),col=c("black","red"),bty="n")
dev.off()

p.VG=sapply(strike, function(x) VG.TimeChange.Euro(param,r,T,sim,S0,x,-1)$price[1]) 
p.GBM=sapply(strike, function(x) GBM.Euro(r,param[1],T,sim,S0,x,-1)$price[1])
png(file="Strike price vs put option price.png",width=7,height=5,units="in",res=300)
plot(strike,p.VG,type="l",xlab="K",ylab="Price",main="Strike price vs. Put option price")
lines(strike,p.GBM,lty=2,col="red")
legend("topleft",legend=c("VG","GBM"),lty=c(1,2),col=c("black","red"),bty="n")
dev.off()


#4 Fit VG distribution
############################################################################################3
## Please see VG_FitDistribution.R


#5 VG distribution vs. different parameters
#########################################################################################
png(file="VG density vs. sigma.png",width=7,height=5,units="in",res=300)
curve(dVG(x,0, theta=0, nu=0.5, sigma=0.25),from=-1,to=1,ylim=c(0,2.5),ylab="density",
      main="VG density vs. sigma")
curve(dVG(x,0, theta=0, nu=0.5, sigma=0.3),add=TRUE,col="blue",lty=2)
curve(dVG(x,0, theta=0, nu=0.5, sigma=0.2),add=TRUE,col="red",lty=3)
legend("topright",legend=c("theta=0,nu=0.5,sigma=0.2",
                           "theta=0,nu=0.5,sigma=0.25",
                           "theta=0,nu=0.5,sigma=0.3"),
       col=c("red","black","blue"),lty=c(3,1,2),bty="n",cex=0.9)
dev.off()

png(file="VG density vs. theta.png",width=7,height=5,units="in",res=300)
curve(dVG(x,0, theta=0, nu=0.5, sigma=0.25),from=-1,to=1,ylab="density",
      main="VG density vs. theta")
curve(dVG(x,0, theta=0.2, nu=0.5, sigma=0.25),add=TRUE,col="blue",lty=2)
curve(dVG(x,0, theta=-0.3, nu=0.5, sigma=0.25),add=TRUE,col="red",lty=3)
legend("topright",legend=c("theta=-0.3,nu=0.5,sigma=0.25",
                           "theta=0,nu=0.5,sigma=0.25",
                           "theta=0.2,nu=0.5,sigma=0.25"),
       col=c("red","black","blue"),lty=c(3,1,2),bty="n",cex=0.9)
dev.off()

png(file="VG density vs. nu.png",width=7,height=5,units="in",res=300)
curve(dVG(x,0, theta=0, nu=0.5, sigma=0.25),from=-1,to=1,ylim=c(0,2.5),ylab="density",
      main="VG density vs. nu")
curve(dVG(x,0, theta=0, nu=0.2, sigma=0.25),add=TRUE,col="red",lty=2)
curve(dVG(x,0, theta=0, nu=0.8, sigma=0.25),add=TRUE,col="blue",lty=3)
legend("topright",legend=c("theta=0,nu=0.2,sigma=0.25",
                           "theta=0,nu=0.5,sigma=0.25",
                           "theta=0,nu=0.8,sigma=0.25"),
       col=c("red","black","blue"),lty=c(2,1,3),bty="n",cex=0.9)
dev.off()

png(file="VG density vs. Time.png",width=7,height=5,units="in",res=300)
curve(dVG(x,0, theta=0, nu=0.5, sigma=0.25,t=1),from=-1,to=1,ylab="density",
      main="VG density vs. Time")
curve(dVG(x,0, theta=0, nu=0.5, sigma=0.25,t=2),add=TRUE,col="blue",lty=2)
curve(dVG(x,0, theta=0, nu=0.5, sigma=0.25,t=3),add=TRUE,col="red",lty=3)
legend("topright",legend=c("t=1",
                           "t=2",
                           "t=3"),
       col=c("black","blue","red"),lty=c(1,2,3),bty="n")
dev.off()


#6 Greeks
################################################################################
delta.diff(param,r,T,sim,S0,K,1,1)
delta.pathwise(param,r,T,sim,S0,K,1)[[1]]
fourier.greeks("delta", S0, K, r, T, param[2] , param[3] , param[1])

gamma.diff(param,r,T,sim,S0,K,1,1)
gamma.pathwise(param,r,T,sim,S0,K,1,1)
fourier.greeks("gamma", S0, K, r, T, param[2] , param[3] , param[1])

rho.diff(param,r,T,sim,S0,K,1,0.1)
rho.pathwise(param,r,T,sim,S0,K,1)
fourier.greeks("rho", S0, K, r, T, param[2] , param[3] , param[1])



sigma.seq=seq(0.01,0.99,by=0.02)
delta=sapply(sigma.seq,function(x) delta.pathwise(c(x,param[2:3]),r,T,sim,S0,K,1)[[1]][1])
gamma=sapply(sigma.seq,function(x) gamma.pathwise(c(x,param[2:3]),r,T,sim,S0,K,1,1)[[1]][1])
rho=sapply(sigma.seq,function(x) rho.pathwise(c(x,param[2:3]),r,T,sim,S0,K,1)[[1]][1])
theta=sapply(sigma.seq,function(x) theta.pathwise(c(x,param[2:3]),r,T,sim,S0,K,1)[[1]][1])
vega=sapply(sigma.seq,function(x) vega.pathwise(c(x,param[2:3]),r,T,sim,S0,K,1)[[1]][1])

delta2=sapply(sigma.seq,function(x) bs.delta(S0, K, r, q = 0 , T, x,1))
gamma2=sapply(sigma.seq,function(x) bs.gamma(S0, K, r, q = 0 , T, x,type = 1))
rho2=sapply(sigma.seq,function(x) bs.rho(S0, K, r, q = 0, T, x, type=1))
vega2=sapply(sigma.seq,function(x) bs.vega(S0, K, r, q = 0, T, x,type = 1))

png(file="Greeks vs. sigma.png",width=7,height=7,units="in",res=300)
par(mfrow=c(2,2))

scatter.smooth(sigma.seq,delta,span=0.2,xlab="sigma",lpars=list(col="red",lwd=2))
lines(sigma.seq,delta2,xlab="sigma",lwd=2,lty=2)
legend("bottomright",c("VG","BS"),col=c("red","Black"),lty=c(1,2),bty="n")

scatter.smooth(sigma.seq,gamma,span=0.2,xlab="sigma",lpars=list(col="red",lwd=2))
lines(sigma.seq,gamma2,xlab="sigma",lwd=2,lty=2)
legend("topright",c("VG","BS"),col=c("red","Black"),lty=c(1,2),bty="n")

scatter.smooth(sigma.seq,rho,span=0.2,xlab="sigma",lpars=list(col="red",lwd=2))
lines(sigma.seq,rho2,xlab="sigma",lwd=2,lty=2)
legend("bottomright",c("VG","BS"),col=c("red","Black"),lty=c(1,2),bty="n")

scatter.smooth(sigma.seq,vega,span=0.2,xlab="sigma",lpars=list(col="red",lwd=2))
lines(sigma.seq,vega2,xlab="sigma",lwd=2,lty=2)
legend("bottomright",c("VG","BS"),col=c("red","Black"),lty=c(1,2),bty="n")

dev.off()


theta.seq=seq(-0.2,0.2,by=0.01)
delta=sapply(theta.seq,function(x) delta.pathwise(c(param[1],x,param[3]),r,T,sim,S0,K,1)[[1]][1])
gamma=sapply(theta.seq,function(x) gamma.pathwise(c(param[1],x,param[3]),r,T,sim,S0,K,1,1)[[1]][1])
rho=sapply(theta.seq,function(x) rho.pathwise(c(param[1],x,param[3]),r,T,sim,S0,K,1)[[1]][1])
theta=sapply(theta.seq,function(x) theta.pathwise(c(param[1],x,param[3]),r,T,sim,S0,K,1)[[1]][1])
vega=sapply(theta.seq,function(x) vega.pathwise(c(param[1],x,param[3]),r,T,sim,S0,K,1)[[1]][1])

png(file="Greeks vs. theta.png",width=7,height=7,units="in",res=300)
par(mfrow=c(2,2))
scatter.smooth(theta.seq,delta,span=0.5,xlab="theta",lpars=list(col="red",lwd=2))
scatter.smooth(theta.seq,gamma,span=0.5,xlab="theta",lpars=list(col="red",lwd=2))
scatter.smooth(theta.seq,rho,span=0.5,xlab="theta",lpars=list(col="red",lwd=2))
#scatter.smooth(theta.seq,theta,span=0.5,xlab="theta",lpars=list(col="red",lwd=2))
scatter.smooth(theta.seq,vega,span=0.5,xlab="theta",lpars=list(col="red",lwd=2))
par(mfrow=c(1,1))
dev.off()


nu.seq=seq(0,1.5,by=0.05)
delta=sapply(nu.seq,function(x) delta.pathwise(c(param[1:2],x),r,T,sim,S0,K,1)[[1]][1])
gamma=sapply(nu.seq,function(x) gamma.pathwise(c(param[1:2],x),r,T,sim,S0,K,1,1)[[1]][1])
rho=sapply(nu.seq,function(x) rho.pathwise(c(param[1:2],x),r,T,sim,S0,K,1)[[1]][1])
theta=sapply(nu.seq,function(x) theta.pathwise(c(param[1:2],x),r,T,sim,S0,K,1)[[1]][1])
vega=sapply(nu.seq,function(x) vega.pathwise(c(param[1:2],x),r,T,sim,S0,K,1)[[1]][1])

png(file="Greeks vs. nu.png",width=7,height=5,units="in",res=300)
par(mfrow=c(2,2))
scatter.smooth(nu.seq,delta,span=0.6,xlab="nu",lpars=list(col="red",lwd=2))
scatter.smooth(nu.seq,gamma,span=0.6,xlab="nu",lpars=list(col="red",lwd=2))
scatter.smooth(nu.seq,rho,span=0.6,xlab="nu",lpars=list(col="red",lwd=2))
#scatter.smooth(nu.seq,theta,span=0.6,xlab="nu",lpars=list(col="red",lwd=2))
scatter.smooth(nu.seq,vega,span=0.6,xlab="nu",lpars=list(col="red",lwd=2))
par(mfrow=c(1,1))
dev.off()


#7 Variance reduction
###########################################################################

## Importance Sampling
param=c(0.25,-0.05,0.5) #c(sigma,theta,nu)
S0=100
K=45
T=1
sim=1e5
r=0.05
type=-1

param.old=param #c(sigma,theta,nu)
param.new=c(0.35,-0.1,0.5)

IS(param.old, param.new,T,type)
VG.TimeChange.Euro(param,r,T,sim,S0,K,type)$price
VG.DiffGamma.Euro(param,r,T,sim,S0,K,type)$price

## Control Variates
control.price(S0, K, r, T, param[1], steps, type=1, param[3], param[2], sim)
vg.geom.asian.mc.price(S0, K, r, T, param[1], type=1, param[3], param[2], sim, steps, "discrete-geometric") 
FFT.price.geom.asian(char.VG, S0 = S0, K = K, r = r, T = T,type = 1, 1000)
