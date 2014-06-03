# Price of a Geometric Asian Call using Analytical Pricing
geom.asian.call <- function(S,K,r,T,vol,N){
  dT  <- T/N
  v   <- r - vol^2/2
  
  a   <- log(S) + v*dT + 0.5*v*(T-dT)
  c   <- vol^2*dT + vol^2*(T-dT)*(2*N-1)/(6*N)
  x   <- (a - log(K) + c) / sqrt(c)
  
  exp(-r*T)*(exp(a+c/2)*pnorm(x) - K*pnorm(x-sqrt(c)))
}

# Compute the stock price evolution for for N steps and M paths
sprice <- function(S, r, T, vol, N, M) {
  tau <- T/N
  mat <- matrix(rnorm( M *  N ), nrow = M)
  val <- exp((r - 0.5*vol^2)*tau + vol*sqrt(tau)*mat)
  val <- cbind(rep(S,M),val)
  t(apply(val,1,cumprod))
}

vg.sim.price <- function(S0, K, r, T, sigma, nu, theta, nsim, steps) {
  tau <- T/steps
  S   <- matrix(S0, nrow = nsim,ncol = 1+steps)
  omega <- 1/nu * log(1 - theta*nu - 0.5*sigma^2*nu)
      
  for (i in 1:steps) {
    t <- rgamma(nsim, shape = tau/nu, scale = nu)
    Z <- rnorm(nsim, 0, 1)
        
    X <- theta*t + sigma*Z*sqrt(t)
    S[,i+1] <- S[,i] * exp(r*tau + omega*tau + X)
  }
  
  S
}

control.price <- function(S, K, r, T, vol, N, type, nu, theta, nsim) {
  
  #pilot simulation
  pilot.size <- round(0.2 * nsim)
  vg.paths   <- vg.sim.price (S0, K, r, T, vol, nu, theta, pilot.size, N)
  bs.paths <- sprice(S,r,T,vol,N,pilot.size)
  
  geom.mean.price  <- exp(apply(log(vg.paths),1,mean))
  log.mean.price  <-  0.5/N*apply(  log(bs.paths[,1:N]) +  
                                      log(bs.paths[,2:(N+1)]),
                                    1, sum )
  
  
  geom.bs.asian.payoff <-  exp(-r*T)*pmax(exp(log.mean.price)- K,0)
  geom.vg.asian.payoff <-  exp(-r*T)*pmax(type*(geom.mean.price- K),0)

  beta    <-  -coef(lm(geom.vg.asian.payoff ~ geom.bs.asian.payoff))
  
  #pricing simulation
  pricing.size <- nsim - pilot.size
  vg.paths   <- vg.sim.price (S0, K, r, T, vol, nu, theta, pricing.size, N)
  bs.paths <- sprice(S,r,T,vol,N,pricing.size)
  
  geom.mean.price  <- exp(apply(log(vg.paths),1,mean))
  log.mean.price  <-  0.5/N*apply(  log(bs.paths[,1:N]) +  
                                      log(bs.paths[,2:(N+1)]),
                                    1, sum )
  
  
  geom.bs.asian.payoff <-  exp(-r*T)*pmax(exp(log.mean.price)- K,0)
  geom.vg.asian.payoff <-  exp(-r*T)*pmax(type*(geom.mean.price- K),0)
  
  geom.asian.price <- geom.asian.call(S,K,r,T,vol,N)
  
  cv.est    <- geom.vg.asian.payoff + beta[2]*(geom.bs.asian.payoff - geom.asian.price)
  cv.price  <- mean(cv.est)
  cv.se     <- sd(cv.est)/pricing.size
  
  list(price = cv.price, se = cv.se)  
}

S0 <- 100; K <- 110; r <- 0.05; T <- 0.25; sigma <- 0.25
theta <- -0.15; nu <- 0.5; nsim <- 10^5; steps <- 10^3

set.seed(882258275)
control.price(S0, K, r, T, sigma, steps, type=1, nu, theta, nsim)

set.seed(403) 
vg.geom.asian.mc.price(S0, K, r, T, sigma, type=1, nu, theta, nsim, steps, "discrete-geometric") 

FFT.price.geom.asian(char.VG, S0 = S0, K = K, r = r, T = T,type = 1, steps)

