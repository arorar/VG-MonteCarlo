# Discretely monitored Fixed Strike Geometric Asian
vg.geom.asian.mc.price <- 
  function(S0, K, r, T, sigma, type, nu, theta, nsim, steps) {

  tau <- T/steps
  S   <- matrix(S0, nrow = nsim,ncol = 1+steps)
  omega <- 1/nu * log(1 - theta*nu - 0.5*sigma^2*nu)
  
  for (i in 1:steps) {
    t <- rgamma(nsim, shape = tau/nu, scale = nu)
    Z <- rnorm(nsim, 0, 1)
  
    X <- theta*t + sigma*Z*sqrt(t)
    S[,i+1] <- S[,i] * exp(r*tau + omega*tau + X)
  }
  
  payoff <- exp(apply(log(S),1,mean)) - K
  prices <- exp(-r*T)*pmax(type*payoff,0) 
  
  list( price = mean(prices), se = sd(prices)/sqrt(nsim))
}

bs.geom.asian.mc.price <- 
  function(S0, K, r, T, sigma, type, nsim, steps) {
    
    tau <- T/steps
    S   <- matrix(S0, nrow = nsim,ncol = 1+steps)
    
    for (i in 1:steps) {
      Z <- rnorm(nsim, 0, 1)
      X <- (r - 0.5*sigma^2)*tau + sigma*sqrt(tau)*Z
      S[,i+1] <- S[,i] * exp(X)
    }
    
    payoff <- exp(apply(log(S),1,mean)) - K
    prices <- exp(-r*T)*pmax(type*payoff,0) 
    
    list( price = mean(prices), se = sd(prices)/sqrt(nsim))
}

S0 <- 100; K <- 110; r <- 0.05; T <- 0.25; sigma <- 0.25
theta <- -0.15; nu <- 0.5; steps <- 100; nsim <- 10^6

mc.put <- bs.geom.asian.mc.price(S0, K, r, T, sigma, type = -1, nsim, steps)
#mc.put <- vg.geom.asian.mc.price(S0, K, r, T, sigma, type = -1, nu, theta, nsim, steps)
