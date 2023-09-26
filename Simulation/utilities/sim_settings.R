# Simulation Settings
# Author : So Young Park

seed <- 1330 # for my simulation
seed0 = 53256

#nsims = 2500
nsims = 1
nsim0=1e5  # for exactLRT() simulation

ss <- seq(0, 1, length.out = 101)
tt <- seq(0, 1, length.out = 41)

# functions for data generation
# fn.phi1 <- function(s){return(rep(1, length(s)))}
# fn.phi2 <- function(s){return(sqrt(2)*sin(2*pi*s))}
fn.phi1 <- function(s){return(sqrt(2)*sin(2*pi*s))}
fn.phi2 <- function(s){return(sqrt(2)*cos(2*pi*s))}

# lam1 <- c(3, 1.5)
# lam2 <- c(2, 1)

lam1 <- c(4,2)
lam2 <- c(3,1)

if(gen.mod == "fpca"){
  fn.xi1.true <- function(t){return( rnorm(1, 0, sqrt(lam1[1]))*sqrt(2)*sin(2*pi*t) + rnorm(1, 0, sqrt(lam1[2]))*sqrt(2)*cos(2*pi*t) )}
  cov.xi1.true <- function(t){
    mat <- lam1[1]*sapply(t, function(a) 2*sin(2*pi*t)*sin(2*pi*a)) +
      lam1[2]*sapply(t, function(a) 2*cos(2*pi*t)*cos(2*pi*a))
    return( mat )
  }

  fn.xi2.true <- function(t){return( rnorm(1, 0, sqrt(lam2[1]))*sqrt(2)*sin(2*pi*t) + rnorm(1, 0, sqrt(lam2[2]))*sqrt(2)*sin(2*pi*t)  )}
  cov.xi2.true <- function(t){
    mat <- lam2[1]*sapply(t, function(a) 2*sin(2*pi*t)*sin(2*pi*a)) + lam2[2]*sapply(t, function(a) 2*sin(2*pi*t)*sin(2*pi*a))
    return( mat )
  }
}

# signal = var(b0) + var(b1)/3 + cov(b0, b1)  = 6, 4
if(gen.mod == "lme"){
  require(mvtnorm)
  beta.vec1.Sigma <- matrix(c(3, 2, 2, 3), nrow = 2)
  beta.vec2.Sigma <- matrix(c(2, 1, 1, 3), nrow = 2)
  if(beta.vec1.Sigma[1,1] + beta.vec1.Sigma[2,2]/3 + beta.vec1.Sigma[1,2] != sum(lam1)){
    warning("recheck sigma (1)")
  }
  if(beta.vec2.Sigma[1,1] + beta.vec2.Sigma[2,2]/3 + beta.vec2.Sigma[1,2] != sum(lam2)){
    warning("recheck sigma (2)")
  }
  
  fn.xi1.true <- function(t){
    beta.vec1 <- rmvnorm(1, mean = rep(0,2), sigma = beta.vec1.Sigma)
    return(beta.vec1 %*% t(cbind(rep(1, length(t)), t)))
    }
  
  fn.xi2.true <- function(t){
    beta.vec2 <- rmvnorm(1, mean = rep(0,2), sigma = beta.vec2.Sigma)
    return(beta.vec2 %*% t(cbind(rep(1, length(t)), t)))
    }
  
  cov.xi1.true <- function(t){
    mat <- cbind(rep(1, length(t)), t) %*% beta.vec1.Sigma %*% t(cbind(rep(1, length(t)), t))
    return( mat )
  }
  cov.xi2.true <- function(t){
    mat <- cbind(rep(1, length(t)), t) %*% beta.vec2.Sigma %*% t(cbind(rep(1, length(t)), t))
    return( mat )
  }
}


# experiment factors :
# mean function
fn.mu.al  <- function(s,t, delta){return(cos(pi*s/2) + delta*t)}
fn.mu.int <- function(s,t, delta){return(cos(pi*s/2) + 5*delta*(t/4 - s)^3)}
fn.mu.c2  <- function(s,t, delta){
  eta_1     <- t
  dd        <- s-0.5; 
  #p1        <- eta_1*{0.5 + 0.5*s + 0.2*s^2  + 0.1*{dd*(dd>0)}^2}
  p1        <- eta_1*{0.05 + 0.5*dd + 0.1*{dd*(dd>0)}}
  p2        <- 0.5 + 0.25*s + 0.1*s^2 + 0.05*{dd*(dd>0)}^2
  return(4*{delta*p1 + p2})
}

