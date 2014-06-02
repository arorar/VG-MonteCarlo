
make.greeks <- function(type, S0, K, r, T, theta , nu , sigma) {
  
  alpha <- 1.65
  i <- sqrt(as.complex(-1))
  
  psi <- function(u) {
    char.VG <- function(u) {
      omega <- 1/nu * log(1 - theta*nu - 0.5*sigma^2*nu)
      exp(i*u*omega*T)*(1 - i*u*nu*theta + 0.5*(sigma*u)^2*nu)^(-T/nu)
    }
    
    phi.risk.neutral <- function(u)  {
      m <- r*T - log(char.VG(-i)) 
      char.VG(u) * exp(i*u*m)
    }
    
    res <- if (type == "delta")  {
      exp(-r*T) * phi.risk.neutral((u - (alpha + 1)*i))/(S0*(i*u + alpha))
    } else if (type == "gamma") {
      exp(-r*T) * phi.risk.neutral((u - (alpha + 1)*i))/S0^2
    } else if (type == "rho") {
      T*exp(-r*T) * phi.risk.neutral((u - (alpha + 1)*i))/(i*u + alpha + 1)
    } else stop ("Other Grreks not available for Levy models")
    
    res
  }
  
  return(psi)
}

fourier.greeks <- function(name, S0, K, r, T, theta , nu , sigma) {
 
  i <- sqrt(as.complex(-1))
  alpha  <- 1.65; N  <- 2^16; h <- 0.05
  
  lambda <- (2 * pi)/(N * h)
  b <- 1/2 * N * lambda
  j <- -b + lambda * (0:(N - 1))
  v <- h * (0:(N - 1))
  
  psi <- make.greeks(name, S0, K, r, T, theta , nu , sigma)
  res <- exp(-alpha * j)/pi * fft( exp(i*b*v) * psi(v) * h/3 * 
                                     (3 + (-1)^(1:N) - as.numeric((1:N) == 1)) )
  inter <- spline(j, Re(res), xout = log(K/S0))
  return(inter$y*S0)
}


##############################
# Finite central difference ##
##############################

delta.diff=function(param,r,T,sim,S0,K,type,h){
  payoff1=VG.TimeChange.Euro(param,r,T,sim,S0+h,K,type)$payoff
  payoff2=VG.TimeChange.Euro(param,r,T,sim,S0-h,K,type)$payoff
  deriv=(payoff1-payoff2)/(2*h)
  deriv.mean=mean(deriv)
  deriv.se=sd(deriv)/sqrt(sim)
  return(c(deriv.mean,deriv.se))
}
#######################################################

gamma.diff=function(param,r,T,sim,S0,K,type,h){
  payoff0=VG.TimeChange.Euro(param,r,T,sim,S0,K,type)$payoff
  payoff1=VG.TimeChange.Euro(param,r,T,sim,S0+h,K,type)$payoff
  payoff2=VG.TimeChange.Euro(param,r,T,sim,S0-h,K,type)$payoff
  deriv=(payoff1+payoff2-2*payoff0)/h^2
  deriv.mean=mean(deriv)
  deriv.se=sd(deriv)/sqrt(sim)
  return(c(deriv.mean,deriv.se))
}
#######################################################

rho.diff=function(param,r,T,sim,S0,K,type,h){
  payoff1=VG.TimeChange.Euro(param,r+h,T,sim,S0,K,type)$payoff
  payoff2=VG.TimeChange.Euro(param,r-h,T,sim,S0,K,type)$payoff
  deriv=((payoff1)-(payoff2))/(2*h)
  deriv.mean=mean(deriv)
  deriv.se=sd(deriv)/sqrt(sim)
  return(c(deriv.mean,deriv.se))
}
####################################################

theta.diff=function(param,r,T,sim,S0,K,type,h){
  param1=c(param[1],param[2]+h,param[3])
  param2=c(param[1],param[2]-h,param[3])
  payoff1=VG.TimeChange.Euro(param1,r,T,sim,S0,K,type)$payoff
  payoff2=VG.TimeChange.Euro(param2,r,T,sim,S0,K,type)$payoff
  deriv=(payoff1-payoff2)/(2*h)
  deriv.mean=mean(deriv)
  deriv.se=sd(deriv)/sqrt(sim)
  return(c(deriv.mean,deriv.se))
}
####################################################

vega.diff=function(param,r,T,sim,S0,K,type,h){
  param1=c(param[1]+h,param[2:3])
  param2=c(param[1]-h,param[2:3])
  payoff1=VG.TimeChange.Euro(param1,r,T,sim,S0,K,type)$payoff
  payoff2=VG.TimeChange.Euro(param2,r,T,sim,S0,K,type)$payoff
  deriv=(payoff1-payoff2)/(2*h)
  deriv.mean=mean(deriv)
  deriv.se=sd(deriv)/sqrt(sim)
  return(c(deriv.mean,deriv.se))
  #return(deriv)
}


################################
# Pathwise #####################
################################

delta.pathwise=function(param,r,T,sim,S0,K,type){
  S=VG.TimeChange.Euro(param,r,T,sim,S0,K,type)$S
  if(type==1){
    v=exp(-r*T)*as.numeric(S>K)*S/S0
  }else{
    v=exp(-r*T)*as.numeric(S<K)*(-S/S0)
  }
  
  deriv.mean=mean(v)
  deriv.se=sd(v)/sqrt(sim)
  return(list(c(deriv.mean,deriv.se),v))
}
#########################################################

gamma.pathwise=function(param,r,T,sim,S0,K,type,h){
  delta1=delta.pathwise(param,r,T,sim,S0+h,K,type)[[2]]
  delta2=delta.pathwise(param,r,T,sim,S0-h,K,type)[[2]]
  v=(delta1-delta2)/(2*h)
  
  deriv.mean=mean(v)
  deriv.se=sd(v)/sqrt(sim)
  return(c(deriv.mean,deriv.se))
}
########################################################

rho.pathwise=function(param,r,T,sim,S0,K,type){
  S=VG.TimeChange.Euro(param,r,T,sim,S0,K,type)$S
  if(type==1){
    v=exp(-r*T)*as.numeric(S>K)*T*K
  }else{
    v=exp(-r*T)*as.numeric(S<K)*(-T*K)
  }
  
  deriv.mean=mean(v)
  deriv.se=sd(v)/sqrt(sim)
  return(c(deriv.mean,deriv.se))
}
#########################################################

theta.pathwise=function(param,r,T,sim,S0,K,type){
  sigma=param[1]
  theta=param[2]
  nu=param[3]
  deriv.w=-1/(1-theta*nu-sigma^2*nu/2)
  output=VG.TimeChange.Euro(param,r,T,sim,S0,K,type)
  S=output$S
  dG=output$dG
  if(type==1){
    v=exp(-r*T)*as.numeric(S>K)*S*(T*deriv.w+dG)
  }else{
    v=exp(-r*T)*as.numeric(S<K)*(-S*(T*deriv.w+dG))
  }
  
  deriv.mean=mean(v)
  deriv.se=sd(v)/sqrt(sim)
  return(c(deriv.mean,deriv.se))
}
###########################################################

vega.pathwise=function(param,r,T,sim,S0,K,type){
  sigma=param[1]
  theta=param[2]
  nu=param[3]
  deriv.w=-sigma/(1-theta*nu-sigma^2*nu/2)
  output=VG.TimeChange.Euro(param,r,T,sim,S0,K,type)
  S=output$S
  dG=output$dG
  Z=output$Z
  if(type==1){
    v=exp(-r*T)*as.numeric(S>K)*S*(T*deriv.w+sqrt(dG)*Z)
  }else{
    v=exp(-r*T)*as.numeric(S<K)*(-S*(T*deriv.w+sqrt(dG)*Z))
  }
  
  deriv.mean=mean(v)
  deriv.se=sd(v)/sqrt(sim)
  return(c(deriv.mean,deriv.se))
}
#################################################################



#Test the result
param=c(0.25,-0.15,0.5)
S0=30; K=30; T=1; sim=1e5; r=0.05

delta.diff(param,r,T,sim,S0,K,1,1)
delta.pathwise(param,r,T,sim,S0,K,1)[[1]]
fourier.greeks("delta", S0, K, r, T, param[2] , param[3] , param[1])

gamma.diff(param,r,T,sim,S0,K,1,1)
gamma.pathwise(param,r,T,sim,S0,K,1,1)
fourier.greeks("gamma", S0, K, r, T, param[2] , param[3] , param[1])

rho.diff(param,r,T,sim,S0,K,1,0.1)
rho.pathwise(param,r,T,sim,S0,K,1)
fourier.greeks("rho", S0, K, r, T, param[2] , param[3] , param[1])

theta.diff(param,r,T,sim,S0,K,1,0.1)
theta.pathwise(param,r,T,sim,S0,K,1)

vega.diff(param,r,T,sim,S0,K,1,0.1)
vega.pathwise(param,r,T,sim,S0,K,1)
