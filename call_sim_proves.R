#Exponencial
print('exponencial')
trans12_exp <- T
scale_param <- 688.1962
shape_param <- 1
pci <- F
fixed_pci <- F
OUT_FILE <- "simul1_exp.RData"
source("sim_proves.R")

#Exponencial as Weibull
print('exponencial as Weibull')
trans12_exp <- F
scale_param <- 688.1962
shape_param <- 1
pci <- F
fixed_pci <- F
OUT_FILE <- "simul1_exp_weib.RData"
source("sim_proves.R")

#Weibull
print('Weibull')
trans12_exp <- F
scale_param <- 198.232490
shape_param <- 1.571168
pci <- T
fixed_pci <- F
OUT_FILE <- "simul1_weibull.RData"
source("sim_proves.R")

#Weibull amb pci fixa
print('Weibull')
trans12_exp <- F
scale_param <- 198.232490
shape_param <- 1.571168
pci <- T
fixed_pci <- T
OUT_FILE <- "simul1_weibull_fix.RData"
source("sim_proves.R")