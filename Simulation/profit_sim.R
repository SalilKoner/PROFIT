

rm(list=ls())

# factor <- 11 ; delta  <- 0.3; PVE <- 0.99;  
# N_Knots <- 20; smooth_mean <- FALSE; ncores <- 4; n_rep <- 20


# IMPORTANT!!!!!! 
args              <- commandArgs(TRUE)
factor            <- as.integer(args[1]) # factor to determine sample size and sparsity level; 
                                         # see "utilities/factor.R", must be an integer between 11 and 21.
smooth_mean       <- as.logical(args[2]) # TRUE for smooth mean, FALSE for non-smooth mean;
                                         # see "utilities/sim_settings.R"., must be logical.
delta             <- as.numeric(args[3]) # The control parameter for departure from null; delta = 0 indicates null.
PVE               <- as.numeric(args[4]) # Percentage of variation explained set to extract the
                                         # eigenfunctions, should be at least 0.9, and less than 1.
N_Knots           <- as.integer(args[5]) # Number of knots used to represent eta_k(t) using truncated
                                         # linear basis, should be at least 10;
n_rep             <- as.integer(args[6]) # Number of replication used to empirically compute the power, 
                                         # should be >= 5000 for delta=0;
ncores            <- as.integer(args[7]) # Number of cores allocated to run the n_rep replications.


for (vars in c("factor", "smooth_mean", "delta", "PVE", "N_Knots", "n_rep", "ncores")){
  cat(paste0(vars, " = ", get(vars), "; "))
}

source("utilities/factor.R")
source("utilities/packages.R")
source("utilities/data_generator.R")
source("utilities/profit.R")
source("utilities/profit_fixed_Fourier.R")

library(conflicted)
conflicts_prefer(dplyr::filter())

# Parallel simulation starts

cl                  <- makeCluster(ncores)
registerDoSNOW(cl)
progress            <- function(nfin, tag) { cat(sprintf('tasks completed: %d; tag: %d\n', nfin, tag)) }
opts                <- list(progress=progress)
packages_req        <- c("refund", "tidyverse", "RLRsim", "nlme", "fda", "rlang")
profit_eigenbasis   <- foreach(i=1:n_rep, .options.snow=opts, .packages =packages_req) %dopar% {
                        dat <- data_generator(delta=delta, n = n, sigma2 = sigma2, m1 = m1, m2 = m2,
                                              sigma2.xi1 = sigma2.xi1, sigma2.xi2 = sigma2.xi2)
                        try(profit_data_driven(Y=dat$Y, visit.time=dat$visit.time, s=dat$s,
                                          m=dat$m, ind.tgrid=dat$ind.tgrid, tt=dat$tt,
                                          pve = PVE, num_knots = N_Knots), 
                            silent = TRUE) } 


# Extracting the minimum of the p-values from obtained from each k, from pLRT.
succ.class.eb           <- sapply(profit_eigenbasis, class) != "try-error"
min_pval.eb             <- profit_eigenbasis[succ.class.eb] %>% map_dbl(~pluck(.x, "min_pValue")[2])
K.eb                    <- profit_eigenbasis[succ.class.eb] %>% map_dbl(~pluck(.x, "K"))
sig.level               <- c(0.01, 0.05, 0.1, 0.15)
power.profit.eigenbasis <- map_dbl(sig.level, ~mean(min_pval.eb < .x/K.eb))
names(power.profit.eigenbasis) <- paste0("\\alpha=", sig.level)


# Creating appropriate folder to save the results
folder_path   <- paste0(sapply(c("smooth_mean", "factor", "delta", "PVE", "N_Knots"), function(i) 
                                   paste0(i, "=",  get(i) )),
                        collapse = .Platform$file.sep)
result_path   <- file.path("Results",  folder_path)
mkdirs(result_path)
cat(result_path, "\n")
save(power.profit.eigenbasis,   file= file.path(result_path,  "power_profit_eigenbasis.Rdata"))



# This part is requested by Reviewer 1: 
# conduct the projection-based test by projecting onto 30 fourier basis. 
# This will take longer time; and it is optional

Fourier <- FALSE  #change Fourier = TRUE if one wants to run the simulation for Fourier
                  # basis as well. Default is FALSE. Only run for eigenbasis. 

if (Fourier){
  packages_req        <- c("refund", "tidyverse", "RLRsim", "nlme", "fda", "rlang")
  profit_fourierbasis <- foreach(i=1:n_rep, .options.snow=opts, .packages =packages_req) %dopar% {
    dat <- data_generator(delta=delta, n = n, sigma2 = sigma2, m1 = m1, m2 = m2,
                          sigma2.xi1 = sigma2.xi1, sigma2.xi2 = sigma2.xi2)
    try(profit_fixed_Fourier(Y=dat$Y, visit.time=dat$visit.time, s=dat$s,
                             m=dat$m, ind.tgrid=dat$ind.tgrid, tt=dat$tt,
                             nb.fourier = 30), 
        silent = TRUE) } 
  succ.class.fb             <- sapply(profit_fourierbasis, class) != "try-error"
  min_pval.fb               <- profit_fourierbasis[succ.class.fb] %>% map_dbl(~pluck(.x, "min_pValue")[2])
  sig.level                 <- c(0.01, 0.05, 0.1, 0.15)
  power.profit.fourierbasis <- map_dbl(sig.level, ~mean(min_pval.fb < .x/31))
  names(power.profit.fourierbasis) <- paste0("\\alpha=", sig.level) 
  save(power.profit.fourierbasis,   file= file.path(result_path,  "power_profit_fourierbasis.Rdata"))
  
}



