rm(list = ls())


vg.sim.price <- function(S0, K, r, T, sigma, type, nu, theta, nsim) {
  t <- rgamma(nsim, shape = T/nu, scale = nu)
  Z <- rnorm(nsim, 0, 1)
  X <- theta*t + Z*sigma*sqrt(t)
  
  omega <- 1/nu * log(1 - theta*nu - 0.5*sigma^2*nu)
  S0 * exp(r*T + omega*T + X)
}

control.price <- function(S0, K, r, T, sigma, type, nu, theta, nsim) {
  
  #pilot simulation
  pilot.size <- round(0.2 * nsim)
  term.price <- vg.sim.price(S0, K, r, T, sigma, type, nu, theta, pilot.size)
  
  fftprice <- sapply(term.price, function(x) FFT.price(char.BS, S0, K, r, T, type))
  bsprice  <- sapply(term.price, function(x) bs.price(x, K, r, T, sigma, type))
  
  mean.price <- S0*exp(r*T)
  disc.payoff <- exp(-r*T)*pmax(type * ( term.price - K) ,0)
  mean.payoff <- mean(disc.payoff)
  
  err.price <-  term.price - mean.price
  err.option.price <- fftprice - bsprice
  err.payoff <- disc.payoff - mean.payoff
  
  fit <- lm(err.payoff ~ -1 + err.price + err.option.price)
  beta <- coef(fit)
  
  #pricing simulation
  pricing.size <- nsim - pilot.size
  term.price <- vg.sim.price(S0, K, r, T, sigma, type, nu, theta, pricing.size)
  
  fftprice <- sapply(term.price, function(x) FFT.price(char.BS, S0, K, r, T, type))
  bsprice  <- sapply(term.price, function(x) bs.price(x, K, r, T, sigma, type))
  
  mean.price <- S0*exp(r*T)
  disc.payoff <- exp(-r*T)*pmax(type * ( term.price - K) ,0)
  mean.payoff <- mean(disc.payoff)
  
  err.price <-  term.price - mean.price
  err.option.price <- fftprice - bsprice
  err.payoff <- disc.payoff - mean.payoff
  
  cv.est    <- disc.payoff - beta[1]*err.price - beta[2]*err.option.price
  cv.price  <- mean(cv.est)
  cv.se     <- sd(cv.est)/pricing.size
  
  list(price = cv.price, se = cv.se)  
}

S0 <- 100; K <- 110; r <- 0.05; T <- 0.25; sigma <- 0.25
theta <- -0.15; nu <- 0.5; nsim <- 10^3

control.price(S0, K, r, T, sigma, type=1, nu, theta, nsim)
