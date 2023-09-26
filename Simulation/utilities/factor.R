
# Moderately sparse case; when m in {8,9,10,11,12}

# n = 100
if(factor == 11){
  gen.mod = "fpca"; source("utilities/sim_settings.R"); source("utilities/meanSwitch.R")
  sigma2 = 10 ; sigma2.xi1 = sum(lam1)/3 ; sigma2.xi2 = sum(lam2)/3
  n = 100; m1 = 8; m2 = 12 
}

# n = 150
if(factor == 12){
  gen.mod = "fpca"; source("utilities/sim_settings.R"); source("utilities/meanSwitch.R")
  sigma2 = 10 ; sigma2.xi1 = sum(lam1)/3 ; sigma2.xi2 = sum(lam2)/3
  n = 150; m1 = 8; m2 = 12
}

# n = 200
if(factor == 13){
  gen.mod = "fpca"; source("utilities/sim_settings.R"); source("utilities/meanSwitch.R")
  sigma2 = 10 ; sigma2.xi1 = sum(lam1)/3 ; sigma2.xi2 = sum(lam2)/3
  n = 200; m1 = 8; m2 = 12
}

# n = 300
if(factor == 14){
  gen.mod = "fpca"; source("utilities/sim_settings.R"); source("utilities/meanSwitch.R")
  sigma2 = 10 ; sigma2.xi1 = sum(lam1)/3 ; sigma2.xi2 = sum(lam2)/3
  n = 300; m1 = 8; m2 = 12
}

# n = 400 (Added later)
if(factor == 19){
  gen.mod = "fpca"; source("utilities/sim_settings.R"); source("utilities/meanSwitch.R")
  sigma2 = 10 ; sigma2.xi1 = sum(lam1)/3 ; sigma2.xi2 = sum(lam2)/3
  n = 400; m1 = 8; m2 = 12
}


# Low sparse case;  when m in {15,16,17,18,19,20}

# n = 100
if(factor == 15){
  gen.mod = "fpca"; source("utilities/sim_settings.R"); source("utilities/meanSwitch.R")
  sigma2 = 10 ; sigma2.xi1 = sum(lam1)/3 ; sigma2.xi2 = sum(lam2)/3
  n = 100; m1 = 15; m2 = 20
}


# n = 150
if(factor == 16){
  gen.mod = "fpca"; source("utilities/sim_settings.R"); source("utilities/meanSwitch.R")
  sigma2 = 10 ; sigma2.xi1 = sum(lam1)/3 ; sigma2.xi2 = sum(lam2)/3
  n = 150; m1 = 15; m2 = 20
}

# n = 200
if(factor == 17){
  gen.mod = "fpca"; source("utilities/sim_settings.R"); source("utilities/meanSwitch.R")
  sigma2 = 10 ; sigma2.xi1 = sum(lam1)/3 ; sigma2.xi2 = sum(lam2)/3
  n = 200; m1 = 15; m2 = 20
}

# n = 300
if(factor == 18){
  gen.mod = "fpca"; source("utilities/sim_settings.R"); source("utilities/meanSwitch.R")
  sigma2 = 10 ; sigma2.xi1 = sum(lam1)/3 ; sigma2.xi2 = sum(lam2)/3
  n = 300; m1 = 15; m2 = 20
}

# n = 400
if(factor == 21){
  gen.mod = "fpca"; source("utilities/sim_settings.R"); source("utilities/meanSwitch.R")
  sigma2 = 10 ; sigma2.xi1 = sum(lam1)/3 ; sigma2.xi2 = sum(lam2)/3
  n = 400; m1 = 15; m2 = 20
}

