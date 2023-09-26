

# factor <- 11 ; delta  <- 0; PVE <- 0.99;  
# N_Knots <- 20; smooth_mean <- TRUE; ncores <- 4; n_rep <- 5000

args              <- commandArgs(TRUE)
factor            <- as.integer(args[1]) # factor to determine sample size and sparsity level; see "utilities/factor.R", must be an integer between 11 and 21.
smooth_mean       <- as.logical(args[2]) # TRUE for smooth mean, FALSE for non-smooth mean; see "utilities/sim_settings.R".=, must be logical.
PVE               <- as.numeric(args[3]) # Percentage of variation explained set to extract the eigenfunctions, should be at least 0.9, and less than 1.
N_Knots           <- as.integer(args[4]) # Number of knots used to represent eta_k(t) using truncated linear basis, should be at least 10;
n_rep             <- as.integer(args[5]) # should be at least 5000 to get verifiable results
ncores            <- as.integer(args[6]) # Number of cores allocated to run the n_rep replications.

delta             <- 0 # Delta must be zero to get the distribution of p-values under null.
if (n_rep < 1000){
  warning("The value of n_rep shoule be atleast 1000 to get reliable QQ plot of pValues")
}

for (vars in c("factor", "smooth_mean", "delta", "PVE", "N_Knots", "n_rep", "ncores")){
  cat(paste0(vars, " = ", get(vars), "; "))
}

source("utilities/factor.R")
source("utilities/packages.R")
source("utilities/data_generator.R")
source("utilities/profit.R")

library(conflicted)
conflicts_prefer(dplyr::filter())

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

succ.class     <- sapply(profit_eigenbasis, class) != "try-error"
pLRT           <- profit_eigenbasis[succ.class]
K.eb           <- profit_eigenbasis[succ.class] %>% map_dbl(~pluck(.x, "K"))
common.K       <- min(K.eb)
pValues        <-  bind_cols(lapply(1:common.K, function(k) pLRT %>% purrr::map_dbl(~pluck(.,"pValue")[2,k]) ) )
names(pValues) <- paste0('Ef', 1:common.K)
pValues_long   <- pValues %>% pivot_longer(everything(), 
                                         names_to = "Ef", values_to = "Pvalue") %>%
                mutate(Ef = factor(Ef))

qq <- pValues_long %>% ggplot(aes(sample=Pvalue, color=Ef)) +
  stat_qq(distribution = qunif, dparams=list(min=0, max=1), size=0.2) +
  stat_qq_line(distribution = qunif, dparams=list(min=0, max=1)) +
  labs(x="Ordered p-values", y="Uniform quantiles") +
  scale_color_manual(values=c("red", "green", "blue", "magenta")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size=1.6), 
        axis.line.x = element_line(size = 1.4),
        axis.text = element_text(face='bold', size=20, color = 'blue'),
        strip.text.x = element_text(face='bold', size=20, color='Red'),
        strip.text.y = element_text(face='bold', size=20, color='Red'),
        strip.background = element_rect(color = "black", fill = 'white', size = 2),
        axis.title = element_text(face = "bold", size=26),
        legend.text=element_text(size=24),
        legend.title=element_text(size=20),
        legend.position = "bottom",
        legend.direction = "horizontal", legend.box = "vertical")

folder_path   <- paste0(sapply(c("smooth_mean", "factor", "delta", "PVE", "N_Knots"), function(i) 
  paste0(i, "=",  get(i) )),
  collapse = .Platform$file.sep)
result_path   <- file.path("Results",  folder_path)
mkdirs(result_path)
cat(result_path, "\n")
ggsave(filename="QQplot_pValues.pdf", plot=qq, path=result_path, 
       device="pdf", width = 11.69, height = 8.27)
