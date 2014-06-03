rm ( list = ls() )

Gamma.bridge <- function(N, rand.g, T, theta, sigma, nu, seed = NULL) {
  
  dt <- T/N; t <- seq(0,T,by=dt)
  X <- rep(NA,N+1); gam <- rep(NA,N+1)
  X[1] <- 0; gam[1] <- 0
  
  gam[N+1] <- rand.g
  X[N+1]   <- rnorm(1,mean=theta*rand.g, sd=sigma*sqrt(rand.g))
  
  for ( k in 1:log(N,2) ) {
    
    n <- 2^(log(N,2) - k)

    for ( j in 1:2^(k - 1) )  {
      
      i <- (2*j - 1)*n + 1
      ind <- which(t == t[i]);
      ind.left <- which(t == t[i-n]);ind.right <- which(t == t[i+n]); 
      
      Y <- rbeta(1, (t[i] - t[i-n])/nu , (t[i+n] - t[i])/nu)
        
      gam[ind] <- gam[ind.left] + (gam[ind.right] - gam[ind.left])*Y
      if (gam[ind.right] <= gam[ind.left]) print(rand.g)
      Z <- rnorm(1, 0, sigma*sqrt((gam[ind.right] - gam[ind])*Y))
      X[ind] <- Y*(X[ind.right]) + (1-Y)*(X[ind.left]) + Z
    }  
  }
  X[-1]
}

asian.payoff <- function(rand.g, S0, K, r, T, vol, steps, type, nu, theta) {
  X <- Gamma.bridge(steps, rand.g, T, theta, vol, nu)
  omega <- 1/nu * log(1 - theta*nu - 0.5*vol^2*nu)
  
  dt <- T/N; S <- rep(NA,N)
  S <- S0 * exp((r + omega)*(1:steps)*dt + X)
  exp(-r*T)*max( type * (exp(mean(log(S))) - K) ,0)
}

stratify <- function(N, n, strata, func, ...) {
  mu.hat <- 0; sigma.sq.hat <- 0
  sigma.strata <- c()
  
  stratas <- length(n)
  
  p <- rep(1/stratas,stratas); 
  X <- qgamma(seq(0,1,length.out=stratas+1),, shape = T/nu, scale = nu)
  
  for (i in 1:stratas) {
    sum <- 0; sum_sq <- 0
    
    if (n[i] == 0) next
    
    for (j in 1:n[i]) {
      
      gam.rnd <- qgamma(pgamma( X[i], shape = T/nu, scale = nu ) + 
                diff( pgamma( X[c(i, i+1)], shape = T/nu, scale = nu ))*runif(1), 
                shape = T/nu, scale = nu)
      
      FUN <- match.fun(func) 
      payoff <- FUN(gam.rnd,...) 
      sum <- sum + payoff
      sum_sq <- sum_sq + payoff^2
    }
    
    mu.hat <- mu.hat + p[i]*sum/n[i]
    sigma.strata[i] <- (sum_sq - sum^2/n[i])/(n[i] - 1)
    sigma.sq.hat <- sigma.sq.hat + ( p[i]* sigma.strata[i] )^2/n[i]
  }
  
  list(price = mu.hat, se = sqrt(sigma.sq.hat), cond.var = sigma.strata)
}

steps <- 64; stratas <- 40; N <- 100
S <- K <- 2; vol <- sigma <- 0.5; T <- 3; r <- 0.05; type <- 1;theta <- -0.15; nu <- 0.25
n <- rep(N,stratas)

pilot <- stratify(N, n, stratas, asian.payoff, S, K, r, T, vol, steps ,type, nu, theta)
p <- rep(1/stratas,stratas); N <- 10^4
n <- round(N*p*sqrt(pilot$cond.var)/sum(p*sqrt(pilot$cond.var)))

print(stratify(N, n, stratas, asian.payoff, S, K, r, T, vol, steps,  type, nu, theta))
