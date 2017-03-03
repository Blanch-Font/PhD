# Simulation 2: Women with constant FP
library(dplyr)
library(msm)
library(survival)
source("aux_fun.R")

# Parameters
# Nsim <- 1000 # number of simulation runs
# N <- 15000 # number of women in the study
# scale_param <- 0.03 # scale parameter for Weibull distribution
aux_scale <- (1/scale_param)^shape_param
# shape_param <- 0.82 # shape parameter for Weibull distribution
# HR_FP <- 2 # theorical hazard ratio
Fup <- 10 # maximum number of visits

msm_50 <- list() # list of MSM models without left-truncation
msm_left <- list() # list of MSM models with left-truncation
cox_50 <- list() # list of Cox models without left-truncation
cox_left <- list() # list of Cox models with left-truncation
disc_50 <- list() # list of Cox models without left-truncation
disc_left <- list() # list of Cox models with left-truncation
msm_50_4 <- list() # list of MSM models without left-truncation
msm_left_4 <- list() # list of MSM models with left-truncation
cox_50_4 <- list() # list of Cox models without left-truncation
cox_left_4 <- list() # list of Cox models with left-truncation
disc_50_4 <- list() # list of Cox models without left-truncation
disc_left_4 <- list() # list of Cox models with left-truncation
for (.nsim in 1:Nsim){
  # Creation of FP as Binomial with Prob(FP)=0.15
  FP <- rbinom(N, 1, 0.15)
  # time to preclinical cancer weibull distribution
  u <- runif(N)
  Tpre <- (-log(u)/(aux_scale*exp(log(HR_FP)*FP)))^(1/shape_param)
  # sojourn time exponential distribution
  u_soj <- runif(N)
  Tsoj <- -log(u_soj)/0.25
  # Clinical time = Time to pre-clinical + sojourn time
  Tclin <- Tpre + Tsoj
  Tenter <- runif(N, 0, 15) # for late entry, Tenter <- 0 if all enter at 50
  # create data set
  # assume fixed screening intervals and every women enters at 50
  my.Data <- data.frame(id = rep(1:N, each = Fup), N_crib = rep(1:Fup, N),
                        FP = rep(FP, each = Fup), Tentry = rep(Tenter, each = Fup),
                        Tpre = rep(Tpre, each = Fup), Tclin = rep(Tclin, each = Fup),
                        Tsoj = rep(Tsoj, each = Fup))
  
  # dataset without left-truncation
  my.Data_50 <- modify_data_2(my.Data)
# datast with left-truncation
  my.Data_left <- modify_data_2(my.Data, T)
  
  ## MSM model
  msm_50[[.nsim]] <- model_msm(my.Data_50)
  msm_left[[.nsim]] <- model_msm(my.Data_left)
  msm_50_4[[.nsim]] <- model_msm(my.Data_50 %>% filter(N_crib <= 4), pci = F)
  msm_left_4[[.nsim]] <- model_msm(my.Data_left %>% filter(N_crib <= 4))
  ## Cox models
  cox_50[[.nsim]] <- model_cox(my.Data_50)
  cox_left[[.nsim]] <- model_cox(my.Data_left)
  cox_50_4[[.nsim]] <- model_cox(my.Data_50 %>% filter(N_crib <= 4))
  cox_left_4[[.nsim]] <- model_cox(my.Data_left %>% filter(N_crib <= 4))
  ## Discret models
  disc_50[[.nsim]] <- model_discret(my.Data_50)
  disc_left[[.nsim]] <- model_discret(my.Data_left)
  disc_50_4[[.nsim]] <- model_discret(my.Data_50 %>% filter(N_crib <= 4))
  disc_left_4[[.nsim]] <- model_discret(my.Data_left %>% filter(N_crib <= 4))
}
save(msm_50, msm_left, cox_50, cox_left, 
     msm_50_4, msm_left_4, cox_50_4, cox_left_4,
     disc_50, disc_left, disc_50_4, disc_left_4,file = OFILE)