bs.price <- function(S0, K, r, T, sigma, type) {
  
  d1 <- (log(S0/K)+(r + sigma^2/2)*T)/(sigma*sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  
  type * ( S0 * pnorm(type * d1) - K*exp(-r*T)*pnorm(type * d2) )
}

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

vg.mc.price <- function(S0, K, r, T, sigma, type, nu, theta, nsim) {
  t <- rgamma(nsim, shape = T/nu, scale = nu)
  Z <- rnorm(nsim, 0, 1)
  X <- theta*t + Z*sigma*sqrt(t)
  
  omega <- 1/nu * log(1 - theta*nu - 0.5*sigma^2*nu)
  S <- S0 * exp(r*T + omega*T + X)
  
  exp(-r*T) * mean(sapply(S, function(x) max(type*(x - K), 0)))  
}

S0 <- 100; K <- 110; r <- 0.05; T <- 0.25; sigma <- 0.25
theta <- -0.15; nu <- 0.5; nsim <- 10^5

FFT.price(char.BS, S0 = S0, K = K, r = r, T = T,type = 1)
FFT.price(char.BS, S0 = S0, K = K, r = r, T = T,type = -1)

bs.price(S0 = S0, K = K, r = r, T = T,sigma = sigma, type = 1)
bs.price(S0 = S0, K = K, r = r, T = T,sigma = sigma, type = -1)

FFT.price(char.VG, S0 = S0, K = K, r = r, T = T, type = 1)
FFT.price(char.VG, S0 = S0, K = K, r = r, T = T, type = -1)

K.seq <- seq(90, 120, length = 50)
mc.price <- c(); fft.price <- c()

for (K in K.seq) {
  mc.put <- vg.mc.price(S0, K, r, T, sigma, type = -1, nu, theta, nsim)
  mc.price <- c(mc.price, mc.put)
  mc.put <- FFT.price(char.VG, S0 , K, r, T, -1)
  fft.price <- c(fft.price, mc.put)
}

png(file = "MC-price-benchmark.png", width = 7, height = 5, units = "in", res = 300)
par(mfrow = c(1,2))
plot(K.seq, 100*abs((mc.price - fft.price)/fft.price), 
     type = "l", xlab = "Strike", ylab = "abs %error")

matplot(K.seq, cbind(mc.price,fft.price), type = "l", 
        xlab = "Strike", col=c("red","blue"),ylab = "price")

legend("topleft", legend=c("mc.price", "fft.price"), 
       col=c("red","blue"), lty=1, bty="n")

par(mfrow = c(1,1))
dev.off()