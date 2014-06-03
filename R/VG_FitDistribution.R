library(moments)
library(quantmod)

dVG <- function(x , c, theta , nu , sigma, t=1) {
  
  temp.func <- function(y) {
    
    integrand <- function(g) {
      dnorm(y-c, mean=theta*g, sd=sigma*sqrt(g))*dgamma(g, shape=t/nu, scale=nu)
    }
    
    integrate(integrand, lower = 0, upper = Inf)$value
  }
  
  sapply(x, temp.func)
}

pVG <- function(x , c, theta , nu , sigma, t=1) {
  
  alpha <- 1
  
  char.VG <- function(u) {
    i <- sqrt(as.complex(-1))
    omega <- 1/nu * log(1 - theta*nu - 0.5*sigma^2*nu)
    exp(i*u*omega*t)*(1 - i*u*nu*theta + 0.5*(sigma*u)^2*nu)^(-t/nu)
  }

  integrand <- function(u) {
    Re(sapply(u, function(y) {
                  i <- sqrt(as.complex(-1))
                  exp(-i*y*x)/(alpha - i*y) * char.VG(y + i*alpha)
              }))
  }
  
  1/(2*pi)*exp(alpha*x)*integrate(integrand, lower=-Inf, upper=Inf,subdivisions=500)$value
}

qVG <- function(p , c, theta , nu , sigma, t=1) {
  
  squared.error <- function(x) {
    (p - pVG(x, c, theta , nu , sigma, t))^2
  }
  
  fit.p <- optim(par=0, fn = squared.error, method = "L-BFGS-B")
  fit.p$par
}

fit.VG <- function(x, trace=FALSE) {
  
  names(x) <- NULL
  x <- coredata(x)
  M <- median(x); S <- skewness(x); K <- kurtosis(x)
  
  sigma <- mad(x)
  nu    <- (K/3 - 1)
  theta <- (S*sigma)/(3*nu)
  c <- M - theta
  
  neg.loglik <- function(x,param) {
    if (trace) print(param)
    -sum(
      log(
        dVG(x, c = param[1], theta = param[2], nu = param[3], sigma = param[4])
      )
    )
  }
  
  start <- c(0, theta, nu, sigma)
  lb <- c(-Inf, -Inf,0,0)
  fit.vg <- optim(par = start, fn = neg.loglik, x = x, method = "L-BFGS-B", lower = lb)
  
  
  list(mean = fit.vg$par[1], theta=fit.vg$par[2], nu=fit.vg$par[3], sigma=fit.vg$par[4])
}

# #download underlying asset prices and calculate its log returns
dat=getSymbols(Symbols="SPY",src="yahoo")
data=Cl(SPY["2012/2014"])
SPreturn=na.omit(diff(log(data)))

#estimate parameters for normal distribution
mu.norm=mean(SPreturn)
sigma.norm=sd(SPreturn)
#estimate robust parameters for normal distribution
mu.norm.rob=median(SPreturn)
sigma.norm.rob=mad(SPreturn)

#estimate parameters for VG distribution
fit=fit.VG(SPreturn)

#compare the density curves of different distributions
png(file="density plots.png",width=7,height=5,units="in",res=300)
plot(density(SPreturn,adjust=0.8),lwd=2,xlab="",
     main="Empirical density vs. Three distributions 
     (sigma=0.0074,theta=-0.001,nu=0.36)")
curve(dnorm(x,mu.norm,sigma.norm),add=TRUE,col="blue",lty=2)
curve(dnorm(x,mu.norm.rob,sigma.norm.rob),add=TRUE,col="green",lty=3)
curve(dVG(x,fit$mean, fit$theta, fit$nu, fit$sigma),add=TRUE,col="red",lwd=2)
legend("topright",legend=c("Empirical density","Normal",
                           "Robust normal","VG"),
       col=c("black","blue","green","red"),lty=c(1,2,3,1),bty="n")
dev.off()


#compare the CDFs
png(file="CDF plots and QQ plot.png",width=9,height=5,units="in",res=300)
par(mfrow=c(1,2))
plot(ecdf(as.numeric(SPreturn)),verticals=TRUE,do.points=FALSE,lwd=2,
     main="Empirical CDF vs. Three distribution CDFs
     (sigma=0.0074,theta=-0.001,nu=0.36)")
curve(pnorm(x,mu.norm,sigma.norm),add=TRUE,col="blue",lty=2)
curve(pnorm(x,mu.norm.rob,sigma.norm.rob),add=TRUE,col="green",lty=3)
curve(pvg(x,vgC=fit$mean, theta=fit$theta, nu=fit$nu, sigma=fit$sigma),add=TRUE,col="red",lwd=2)
legend("topleft",legend=c("Empirical CDF","Normal CDF",
                           "Robust Normal CDF","VG CDF"),
       col=c("black","blue","green","red"),lty=c(1,2,3,1),bty="n",cex=0.8)

qqnorm(SPreturn)
qqline(SPreturn)
dev.off()