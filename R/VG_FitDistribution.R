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
fit=vgFit(SPreturn,startValues="US",paramStart=c(mu.norm.rob,sigma.norm.rob,0,0.5))
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
