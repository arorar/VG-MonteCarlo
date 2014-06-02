#c(sigma,theta,nu)

VG.density<-function(x,sigma,theta,nu){
  v=1/(nu*abs(x))*exp(theta/sigma^2*x-1/sigma*sqrt(2/nu+theta^2/sigma^2)*abs(x))
}

integrand<-function(x,param.old, param.new){

  k.old=VG.density(x,param.old[1],param.old[2],param.old[3])
  k.new=VG.density(x,param.new[1],param.new[2],param.old[3])
  
  k.old-k.new
}

VG.mu<-function(param,sign){
  v=(sqrt(param[2]^2+2*param[1]^2/param[3])+sign*param[2])/2
}

IS <- function(param.old, param.new,T,type) {
  output=VG.DiffGamma.Euro(param.new,r,T,sim,S0,K,type)
  
  mu.old.p=VG.mu(param.old,1)
  mu.old.m=VG.mu(param.old,-1)
  mu.new.p=VG.mu(param.new,1)
  mu.new.m=VG.mu(param.new,-1)
  
  g.p=output$g.p
  g.m=output$g.m
  
  phi.p=exp(2*(mu.new.m/param.new[1]^2-mu.old.m/(param.old[1]^2))*abs(g.p))
  phi.m=exp(2*(mu.new.p/param.new[1]^2-mu.old.p/(param.old[1]^2))*abs(g.m))

#   phi.p=exp(-(mu.new.p*nu-mu.old.p*nu)*abs(g.p))
#   phi.m=exp(-(mu.new.m*nu-mu.old.m*nu)*abs(g.m))
  
  Z=integrate(function(y) integrand(y,param.old, param.new), -Inf,Inf)
  rad.nik=exp(-T * Z$value)* phi.p * phi.m
  
  
  payoff=output$payoff*rad.nik
  p.hat=mean(payoff)
  p.se=sd(payoff)/sqrt(sim)
  return(c(price=p.hat,se=p.se))
}


param.old=c(0.25,0,0.5)
param.new=c(0.4,-0.05,0.5)
IS(param.old, param.new,T,-1)

# sigma=param[1]
# theta=param[2]
# nu=param[3]
# w=log(1-theta*nu-sigma^2*nu/2)/nu
# log(K/S0)-(r+w)*T
