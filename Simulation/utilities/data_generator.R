# data genrating function
#   - generated t_ij by sparsifying the regular grid
#   - sigma2.xi1 == sigma2.xi2
# Author : So Young Park 

data_generator <- function(delta = 0, n = 200, m1 = 10, m2 = 15,
                           sigma2 = 0.5, sigma2.xi1 = 0.5, sigma2.xi2 = 0.5){

  m <- c()
  visit.time <- X <- ind.samp_ <- list()
  for(i in 1:n){
    if(m1 == m2){
      mi <- m1
    }else{
      mi <- sample(seq(m1,m2), size = 1)
    } 
    
    ind.samp <- ind.samp_[[i]]<- sort(sample(1:length(tt), size = mi, replace = FALSE))
    tij <- tt[ind.samp]
    
    xi1 <- fn.xi1.true(tij) + rnorm(mi, 0, sqrt(sigma2.xi1))
    xi2 <- fn.xi2.true(tij) + rnorm(mi, 0, sqrt(sigma2.xi2))   
    XI <- cbind(as.vector(xi1), as.vector(xi2))      


    PHI <- cbind(fn.phi1(ss), fn.phi2(ss))
    
    MU <- t(sapply(tij, function(a) fn.mu(s = ss, t = a, delta = delta )))
    
    visit.time[[i]] <- tij
    m <- c(m, mi)
    X[[i]] <- MU + XI %*% t(PHI)
  }
  
  visit.time <- unlist(visit.time)
  X <- do.call(rbind, X)
  Error <- rnorm(length(X), mean = 0, sd = sqrt(sigma2))
  Y <- X + matrix(Error, nrow=nrow(X))
  
  return(list(m = m, s = ss, visit.time = visit.time, X = X, Y = Y,
              Phi = PHI, ind.tgrid = ind.samp_, tt = tt))
}