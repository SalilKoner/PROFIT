rm(list=ls())

smooth_mean.all        <- list(FALSE, TRUE)
factor.all             <- list(11, 12, 13, 14, 15, 16, 17, 18)
delta.all              <- list(0, 0.05, 0.07, 0.09, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6)
PVE.all                <- list(0.9, 0.95, 0.99)
order_linear_basis.all <- list(1, 2, 3)
N_Knots.all            <- list(10, 20, 30)


arguments_list         <- ls()[stringr::str_detect(ls(), ".all")]
sim.settings           <- lapply(arguments_list, get)
names(sim.settings)    <- stringr::str_remove_all(arguments_list, ".all")
sim.design.all         <- tidyr::expand_grid(!!!sim.settings) 


save(sim.design.all, file = "sim_design.Rdata")


# load("sim_design.Rdata")  
# abc <- cbind(sim.design.all$smooth_mean==FALSE, sim.design.all$delta==0,
#                sim.design.all$factor %in% c(11,13,15,17),
#                sim.design.all$N_Knots==20,
#                sim.design.all$order_linear_basis == 1,
#                sim.design.all$PVE==0.99 )
#   which(apply(abc,1,prod)==1)
  