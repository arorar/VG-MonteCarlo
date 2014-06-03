library(VarianceGamma)

VG.CGM <- function(C,G,M) {
  nu <- 1/C
  theta <- C*(1/M - 1/G)
  sigma <- sqrt(C/2*((1/M + 1/G)^2-(1/M - 1/G)^2))
  c(0, sigma, theta, nu)
}


integrand <- function(x, C, G.old, M.old, G.new, M.new) {
  param <- VG.CGM(C, G.old, M.old)
  d.new <- dvg(x, param = param)
  
  param <- VG.CGM(C, G.new, M.new)
  d.old <- dvg(x,param = param)
  
  d.new - d.old
}

rad.nik <- function(t, C, G.old, M.old, G.new, M.new) {
  g.p <- rgamma(1,shape = C*t, rate = M.old )
  g.m <- rgamma(1,shape = C*t, rate = G.old )

  phi.m <- exp(-(G.new - G.old)*abs(-g.m))
  phi.p <- exp(-(M.new - M.old)*g.p)

  Z <- integrate(integrand, -Inf,Inf, C, G.old, M.old, G.new, M.new)
  exp(-t * Z$value)* phi.p * phi.m
}

 
