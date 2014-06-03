# Discretely monitored Fixed Strike Geometric Asian

asian.discpayoff <- function(id,S, K, r, q, T, vol, N, M, type) {
  price <- sprice(S = S, r = r, q = q, T = T, vol = vol, N = N, M = M)
  avgprice <-  
  exp(-r*T)*pmax( type * (avgprice - K),0)
}


vg.geom.asian.mc.price <- 
  function(S0, K, r, T, sigma, type, nu, theta, nsim, steps, asiantype) {

  tau <- T/steps
  S   <- matrix(S0, nrow = nsim,ncol = 1+steps)
  omega <- 1/nu * log(1 - theta*nu - 0.5*sigma^2*nu)
  
  for (i in 1:steps) {
    t <- rgamma(nsim, shape = tau/nu, scale = nu)
    Z <- rnorm(nsim, 0, 1)
  
    X <- theta*t + sigma*Z*sqrt(t)
    S[,i+1] <- S[,i] * exp(r*tau + omega*tau + X)
  }
  
  payoff <- if(asiantype == "discrete-geometric") {
    exp(apply(log(S),1,mean)) - K
  } else if (asiantype == "discrete-arithmetic") {
    apply(S,1,mean) - K
  }else if (asiantype == "continuous-arithmetic") {
    0.5/steps*apply(S[,1:steps] +  S[,2:(steps+1)],1,sum) - K
  }
  
  prices <- exp(-r*T)*pmax(type*payoff,0) 
  
  list( price = c(price=mean(prices), se = sd(prices)/sqrt(nsim)))
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
    
    list( price = c(price=mean(prices), se = sd(prices)/sqrt(nsim)),S=S)
}

FFT.price.geom.asian <- function(phi, S0, K, r, T, type, steps) {
  
  i <- sqrt(as.complex(-1))
  alpha  <- type*1.65; N  <- 2^16; h <- 0.05
  
  lambda <- (2 * pi)/(N * h)
  b <- 1/2 * N * lambda
  j <- -b + lambda * (0:(N - 1))
  v <- h * (0:(N - 1))
  
  tau <- T/steps
  
  phi.risk.neutral <- function(u)  {
    prod <- 1
    
    for (j in 1:steps){
      m <- r*tau - log(phi(-i)) 
      fact <- (steps - j + 1)/(steps + 1)
      prod <- prod *phi(u*fact)
    }
  
    exp(i*u*m*steps/2) * prod
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
  exp(-0.5 * tau * (sigma*u)^2)
}

char.VG <- function(u) {
  i <- sqrt(as.complex(-1))
  omega <- 1/nu * log(1 - theta*nu - 0.5*sigma^2*nu)
  (1 - i*u*nu*theta + 0.5*(sigma*u)^2*nu)^(-tau/nu)
}


S0 <- 100; K <- 110; r <- 0.05; T <- 0.25; sigma <- 0.25
theta <- -0.15; nu <- 0.5; nsim <- 10^5; steps <- 10^3; tau <- T/steps

FFT.price.geom.asian(char.BS, S0 = S0, K = K, r = r, T = T,type = 1, steps)
mc.call <- bs.geom.asian.mc.price(S0, K, r, T, sigma, type = 1, nsim, steps)
FFT.price.geom.asian(char.VG, S0 = S0, K = K, r = r, T = T,type = 1, steps)
mc.vg.call <- vg.geom.asian.mc.price(S0, K, r, T, sigma, type = 1, nu, theta, nsim, steps)
