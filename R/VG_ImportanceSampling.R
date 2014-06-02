#c(sigma,theta,nu)

VG.density<-function(x,sigma,theta,nu){
  v=1/(nu*abs(x))*exp(theta/sigma^2*x-1/sigma*sqrt(2/nu+theta^2/sigma^2)*abs(x))
}

integrand<-function(x,param.old, param.new){
  nu=param.old[3]
  sigma.old=param.old[1]
  theta.old=param.old[2]
  sigma.new=param.new[1]
  theta.new=param.new[2]
  k.old=VG.density(x,sigma.old,theta.old,nu)
  k.new=VG.density(x,sigma.new,theta.new,nu)
  
  k.old-k.new
}

VG.mu<-function(sigma,theta,nu,sign){
  v=(sqrt(theta^2+2*sigma^2/nu)+sign*theta)/2
}

IS <- function(param.old, param.new,T,type) {
  output=VG.DiffGamma.Euro(param.new,r,T,sim,S0,K,type)
  
  nu=param.old[3]
  sigma.old=param.old[1]
  theta.old=param.old[2]
  sigma.new=param.new[1]
  theta.new=param.new[2]  
  
  mu.old.p=VG.mu(sigma.old,theta.old,nu,1)
  mu.old.m=VG.mu(sigma.old,theta.old,nu,-1)
  mu.new.p=VG.mu(sigma.new,theta.new,nu,1)
  mu.new.m=VG.mu(sigma.new,theta.new,nu,-1)
  
  g.p=output$g.p
  g.m=output$g.m
  
  phi.p=exp(2*(mu.new.m/sigma.new^2-mu.old.m/(sigma.old^2))*abs(g.p))
  phi.m=exp(2*(mu.new.p/sigma.new^2-mu.old.p/(sigma.old^2))*abs(g.m))

#   phi.p=exp(-(mu.new.p*nu-mu.old.p*nu)*abs(g.p))
#   phi.m=exp(-(mu.new.m*nu-mu.old.m*nu)*abs(g.m))
  
  Z=integrate(function(y) integrand(y,param.old, param.new), -Inf,Inf)
  rad.nik=exp(-T * Z$value)* phi.p * phi.m
  
  
  h=output$payoff
  payoff=h*rad.nik
  p.hat=mean(payoff)
  p.se=sd(payoff)/sqrt(sim)
  return(c(price=p.hat,se=p.se))
}sigma=param[1]



# sigma.seq=seq(0,0.5,by=0.01)
# theta.seq=seq(-0.5,0.5,by=0.01)

param.old=c(0.25,0,0.5)
param.new=c(0.4,-0.05,0.5)
IS(param.old, param.new,T,-1)

# sigma=param[1]
# theta=param[2]
# nu=param[3]
# w=log(1-theta*nu-sigma^2*nu/2)/nu
# log(K/S0)-(r+w)*T
