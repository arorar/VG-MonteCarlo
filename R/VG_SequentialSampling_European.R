
#####################################################################
# Price European option via VG Gamma Time-Changed Brownian Motion ###
#####################################################################
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


############################################################
## Price European option via VG Difference of Gammas #######
############################################################
VG.DiffGamma.Euro=function(param,r,T,sim,S0,K,type){
  dt=T
  sigma=param[1]
  theta=param[2]
  nu=param[3]
  w=log(1-theta*nu-sigma^2*nu/2)/nu
  mu.p=VG.mu(param,1)
  mu.m=VG.mu(param,-1)
  nu.p=mu.p^2*nu
  nu.m=mu.m^2*nu
  
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
                        X=X,
                        S=S,
                        g.p=g.p,
                        g.m=g.m))
}

######################################################
## Price European option via GBM Monte Carlo #########
######################################################
GBM.Euro=function(r,sigma,T,sim,S0,K,type){
  dt=T
  eps=rnorm(sim)
  S=S0*exp((r-0.5*sigma^2)*dt+sigma*sqrt(dt)*eps)
  payoff.disc=pmax(type*(S-K),0)*exp(-r*T)
  p.hat=mean(payoff.disc)
  p.se=sd(payoff.disc)/sqrt(sim)
  return(list(price=c(price=p.hat,se=p.se),payoff=payoff.disc,S=S))
}


###################################################################
## Price European option via BS formula ###########################
###################################################################
bs.price <- function(S0, K, r, T, sigma, type) {
  
  d1 <- (log(S0/K)+(r + sigma^2/2)*T)/(sigma*sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  
  type * ( S0 * pnorm(type * d1) - K*exp(-r*T)*pnorm(type * d2) )
}


#####################################################################
## Price European option via Fast Fourier Transform #################
#####################################################################
FFT.price <- function(phi, S0, K, r, T, type) {
  
  i <- sqrt(as.complex(-1))
  alpha  <- type*1.65; N  <- 2^16; h <- 0.05
  
  lambda <- (2 * pi)/(N * h)
  b <- 1/2 * N * lambda
  j <- -b + lambda * (0:(N - 1))
  v <- h * (0:(N - 1))
  
  phi.risk.neutral <- function(u)  {
    m <- r*T - log(phi(-i)) 
    phi(u) * exp(i*u*m)
  }
  
  psi <- function(u) 
    exp(-r*T) * phi.risk.neutral((u - (alpha + 1)*i))/((i*u + alpha)*(i*u + alpha + 1))
  
  res <- exp(-alpha * j)/pi * fft( exp(i*b*v) * psi(v) * h/3 * 
                                     (3 + (-1)^(1:N) - as.numeric((1:N) == 1)) )
  inter <- spline(j, Re(res), xout = log(K/S0))
  return(inter$y*S0)
}

char.BS <- function(u)  { 
  i <- sqrt(as.complex(-1))
  exp(T*(i*u*(-0.5 * sigma^2) - 0.5 * (sigma*u)^2))
}

char.VG <- function(u) {
  i <- sqrt(as.complex(-1))
  omega <- 1/nu * log(1 - theta*nu - 0.5*sigma^2*nu)
  exp(i*u*omega*T)*(1 - i*u*nu*theta + 0.5*(sigma*u)^2*nu)^(-T/nu)
}

#########################################################################################

# #Test the result
# param=c(0.25,-0.05,0.5) #c(sigma,theta,nu)
# S0=100
# K=110
# T=0.25
# sim=1e5
# r=0.05
# type=1
# 
# 
# VG.TimeChange.Euro(param,r,T,sim,S0,K,type)$price
# VG.DiffGamma.Euro(param,r,T,sim,S0,K,type)$price
# GBM.Euro(r,param[1],T,sim,S0,K,type)$price
# bs.price(S0, K, r, T, param[1], type) 

