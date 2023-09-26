
#' PROFIT: Projection based Functional Invariance Testing Algorithm 
#'
#' This algorithm uses the data-driven eigen functions from the marginal covariance function estimated by
#' fast sandwich smoother proposed by Xiao et al (2013). The null distribution of the test statistic is 
#' approximated by RLRSim package.
#' @author Salil Koner \cr Maintainer: Salil Koner
#' \email{skoner@@ncsu.edu}
#' @import refund
#' @import nlme
#' @import RLRsim
#' @param Y : A N by R matrix where each row corresponds to functional trajectories for 
#' subjects across different visit times. N = sum_i m_i
#' @param s : A vector of length R corresponding to the grid points where the trajectories Y_{ij}(s) are observed
#' @param tt : The vector of unique time points $t_{ij}$ where the response Y_{ij}(s) are observed, 
#' pooled over all subjects
#' @param m : A vector of length n where each element corresponds to number of visits times for each subject.
#' this should match sapply(visit.time, length)
#' @param visit.time : A list of length n where each element is a vector of length m[i] indicating the visit times 
#' for each subject. The vectors in visit.time must be a subset of the vector tt.
#' @param ind.tgrid : A list of length n. Each elements is a m[i] length vector indicating the index points 
#' of the visit times. For example, visit.time[[i]] = tt[ind.tgrid[[i]]] for all i=1,2..,n
#' @return The factorization A and a list of weights matrix R with length of the list as the length of X
#' @export
#' @examples
#' factor <- 11; delta <- 0;
#' source("utilities/factor.R")
#' source("utilities/packages.R")
#' source("utilities/data_generator.R")
#' source("utilities/profit.R")
#' require(nlme) ; require(tidyverse) ; require(RLRsim) ; require(refund) ;
#' dat <- data_generator(delta=delta, n = n, sigma2 = sigma2, m1 = m1, m2 = m2,
#'                       sigma2.xi1 = sigma2.xi1, sigma2.xi2 = sigma2.xi2)
#' result <- profit_data_driven(Y=dat$Y, visit.time=dat$visit.time, s=dat$s,
#'                         m=dat$m, ind.tgrid=dat$ind.tgrid, tt=dat$tt)   
#' print(result)                                            


# factor <- 17; smooth_mean <- TRUE; delta <- 0.2; num_knots <- 20; degree_pol <- 1; pve <- 0.99
# filename <- paste0('factor', factor)
# source("utilities/factor.R")
# source("utilities/packages.R")
# source("utilities/data_generator.R")
# dat <- data_generator(delta=delta, n = n, sigma2 = sigma2, m1 = m1, m2 = m2,
#                       sigma2.xi1 = sigma2.xi1, sigma2.xi2 = sigma2.xi2)
# Y=dat$Y; visit.time=dat$visit.time; s=dat$s; m=dat$m; ind.tgrid=dat$ind.tgrid; tt=dat$tt;


profit_data_driven      <- function(Y, s, tt, m, visit.time, ind.tgrid, 
                               pve=0.9, num_knots=20){
  
  nt               <- length(tt)
  R                <- length(s)
  n                <- length(m)
  id               <- unlist(purrr::map(1:n, ~ rep(.x, m[.x])))
  
  not.represent    <- setdiff(1:nt, unique(unlist(ind.tgrid)))
  
  
  # Part 1: Estimation Step for the bivariate mean function and eigen functions of marginal covariance
  
  Ymat.fit         <- fbps(Y, covariates=list(x=visit.time, z=s)) # Estimating the bivariate mean function. 
                                                                  # One can use a different bivariate smoother here
  Ymat.demeaned    <- Y - Ymat.fit$Yhat                           # Computing the residuals
  Ymat.fpca        <- fpca.face(Ymat.demeaned, pve=pve)           # Estimated eigen functions for the marginal FPCA
  Phi.hat          <- (Ymat.fpca$efunctions)*sqrt(R)              # Scaling the eigen functions
  npc              <- Ymat.fpca$npc                               # Projection onto the space of eigen functions. This is a N by npc matrix
  Ymat.proj        <- Y%*%Phi.hat/R
  
  # Part 2: Testing step in the projected model
  p                <- 1                                           # Truncated linear basis p=1
  K                <- num_knots                                          # Number of knot points in truncated linear basis expansion of eta_k(t)
  knots            <- unname(quantile(tt,probs=seq(0,1, length = K+2)[-c(1,K+2)])) # knot points
  #X                <- cbind(1, visit.time)                        # X matrix for linear mixed model
  X                <- do.call(cbind, lapply(0:p, function(i) visit.time^i )) # Added by SK on Aug 23, 2023
  Z                <- outer(visit.time, knots, FUN="-")           # Z matrix for linear mixed model
  Z                <- {Z*(Z > 0)}^p  
  
  TestStat         <- matrix(NA, 2, npc)                      
  pValue           <- matrix(NA, 2, npc)
  
  for (irow in 1:npc){ # Each irow represents the testing corresponding to each of the npc eigen functions

    y                <- Ymat.proj[,irow] # the projected response W_{k,ij}
    dat.fpca         <- data.frame(id=c(rep(0,length(not.represent)), id), 
                                   tindex=c(not.represent, unlist(ind.tgrid)), value=y) %>%
                        spread(tindex, value) %>% dplyr::filter(id > 0) %>% dplyr::select(-id) %>% as.matrix
    
    # Estimation of the covariance matrix Sigma_k using the proxy projections
    fpcaObj          <- fpca.sc(Y=dat.fpca, pve=0.9)  
    
    # Extract the eigen functions for gamma_k(t,t')
    ef               <- fpcaObj$efunctions
    ef.all           <- do.call("rbind", lapply(ind.tgrid, function(t) ef[t,,drop=FALSE]))
    
    dat_lme          <- cbind.data.frame(rep(1,length(id)), id, y, X, Z, ef.all)
    names(dat_lme)   <- c("one", "id", "Y", paste0("X",1:ncol(X)),
                         paste0("Z",1:ncol(Z)), paste0("EF",1:ncol(ef.all)))
    
    random.model     <- as.formula(paste0("~ -1 +",paste0("Z",1:ncol(Z), collapse = "+")))
    ef.model         <- as.formula(paste0("~ -1 +",paste0("EF",1:ncol(ef.all), collapse = "+")))
    
    # Added by SK on Aug 23, 2023
    fixed.model      <- as.formula(paste0("Y ~ -1 + ", paste0("X",1:ncol(X), collapse = " + ")))
    #print(fixed.model)
    lmeObj           <- lme(fixed.model,
                            random = list(one = pdIdent(random.model),
                                          id  = ef.model),
                            data = dat_lme , control  = lmeControl(opt='optim')
    )
    
    # lmeObj           <- lme(Y ~ -1 + X1 + X2,
    #                        random = list(one = pdIdent(random.model),
    #                                      id  = ef.model),
    #                        data = dat_lme , control  = lmeControl(opt='optim')
    # )
    
    # Estimate the eigen values and the nugget using a linear mixed model approach
    sig             <- lmeObj$sigma^2
    m               <- lapply(rev(lmeObj$modelStruct$reStruct), VarCorr)
    ev              <- as.vector(m$id[,"Variance"])*sig
    
    
    # Estimaion of Covariance matrix Sigma_ik for each each subject
    vcov            <- tcrossprod(ef * rep(sqrt(ev), rep.int(nrow(ef), length(ev))))
    diag(vcov)      <- diag(vcov) + ifelse(sig > 0, sig, 1e-5)
    sighat          <- Matrix::bdiag(purrr::map(ind.tgrid, ~ `[`(vcov, i=.x,j=.x))) # Estimated block diagonal matrix \widehat{\Sigma}_k
    
    
    chol_sig        <- t(chol(sighat))
    X.hat           <- forwardsolve(chol_sig, X) # Pre-whiten the Design matrix X
    Z.hat           <- forwardsolve(chol_sig, Z) # Pre-whiten the Design matrix Z
    Y.hat           <- as.vector(forwardsolve(chol_sig, y)) # Pre-whiten the response y
    
    # construction of pseudo-likelihood test-statistic
    
    
    # Method 1: Using NLME package to set up the alternate model. Then compute the test-statistic along with
    # generating sample from the null distribution using RLRSim package
    modeldat        <- cbind.data.frame(Y.hat, X.hat, Z.hat)
    names(modeldat) <- c("Y", paste0("X",1:ncol(X.hat)), paste0("Z",1:ncol(Z.hat)))
    modeldat$one    <- rep(1, nrow(modeldat))
    random.model    <- paste0("~ -1 +",paste0("Z",1:ncol(Z.hat), collapse = "+"))
    
    # Added by SK on Aug 23, 2023
    mA_model <- rlang::call2("lme", !!!list(fixed=reformulate(termlabels = paste0("X",1:ncol(X.hat)), 
                                                              response = 'Y', intercept=F),
                                random=list(one=pdIdent(as.formula(random.model))),
                                data=quote(modeldat), method="ML"))
    mA <- eval(mA_model)
    # mA              <- lme(Y ~ X1 + X2 - 1, random=list(one=pdIdent(as.formula(random.model))),
    #                        data=modeldat, method="ML") # Alternate model
    m0              <- lm(Y ~ X1 - 1, data=modeldat) # null model
    eLRT            <- RLRsim::exactLRT(m=mA,m0=m0,nsim=1e5) # pvalue and the test-statistic
    
    # Method 2: Using direct optimzation to estimate the components of the null
    # and alternate model and compute the test-statistic. Then use RLRSim just to generate s
    # sample from null distribution
    XhtXh       <- crossprod(X.hat)
    ZhtZh       <- crossprod(Z.hat)
    ZhtXh       <- crossprod(Z.hat,X.hat)
    XhtYh       <- as.vector(crossprod(X.hat,Y.hat))
    ZhtYh       <- as.vector(crossprod(Z.hat,Y.hat))
    YhtYh       <- sum(Y.hat^2)
    eig.ZhtZh   <- svd(ZhtZh)$d
    
    
    neg_profile_log_lik <- function(loglam, XtX, ZtX, XtY, ZtZ, ZtY, YtY, eig.ZtZ){
      lam          <- exp(loglam)
      H            <- lam*ZtZ
      diag(H)      <- diag(H) + 1
      chol.H       <- chol(H)
      choltinvZtX  <- forwardsolve(t(chol.H), ZtX)
      XtHinvX      <- XtX - lam*crossprod(choltinvZtX)
      choltinvZtY  <- as.vector(forwardsolve(t(chol.H), ZtY))
      XtHinvY      <- XtY - lam*as.vector(crossprod(choltinvZtX, choltinvZtY))
      chol.XtHinvX <- chol(XtHinvX)
      p2.lik       <- as.vector(forwardsolve(t(chol.XtHinvX), XtHinvY))
      p2.lik2      <- sum(p2.lik^2)
      p1.lik2      <- YtY - lam*sum(choltinvZtY^2)
      p1.lik2 - p2.lik2 + sum(log(1+lam*eig.ZtZ))
    }
    
    neglikfit      <- optimize(f=neg_profile_log_lik, interval=c(-40,40),
                               XtX=XhtXh, ZtX=ZhtXh, XtY=XhtYh, ZtZ=ZhtZh,
                               ZtY=ZhtYh, YtY=YhtYh, eig.ZtZ=eig.ZhtZh)
    
    negLik.H0uHA   <- neglikfit$objective
    negLik.H0      <- YhtYh - (1/XhtXh[1,1])*((XhtYh[1])^2)
    pLRT           <- negLik.H0 - negLik.H0uHA
    sim.LRT        <- RLRsim::LRTSim(X=X.hat, Z=Z.hat, q=ncol(X.hat)-1, sqrt.Sigma=diag(ncol(Z.hat)))
    pval.LRT       <- mean(sim.LRT > pLRT)
    
    pValue[, irow] <- c(eLRT$p.value, pval.LRT) # Both the p-values eLRT$p.value and pval.LRT should be close.
    TestStat[,irow]<- c(eLRT$statistic,pLRT)    # Both the test-statistic eLRT$statistic (computed through RLRSim) and 
                                                # pLRT (obtained by direct optimization) should be close
    
  }
  
  list("max_TestStat"    = apply(round(TestStat,3),1,max), 
       "pValue"          = round(pValue,3),
       "min_pValue"      = apply(round(pValue,3),1,min),
       "K"               = npc,
       "K_with_min_pVal" = apply(round(pValue,3),1,which.min) 
  )
  
}







