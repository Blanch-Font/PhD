# Simulation 4: Women with a time-depending FP as step function and random censoring
library(dplyr)
library(survival)
library(msm)

source("aux_fun.R")

# Parameters
aux_scale <- (1/scale_param)^shape_param
Fup <- 10 # maximum number of visits

msm_50 <- list() # list of MSM models without left-truncation
for (.nsim in 1:Nsim){
  # set.seed(n_seed[.nsim])
  # Uniform for weibull distribution
  u <- runif(N)
  # sojourn time exponential distribution
  u_soj <- runif(N)
  Tsoj <- -log(u_soj)/0.25
  Tenter <- runif(N, 0, 15) # for late entry, Tenter <- 0 if all enter at 50
  # create data set
  # assume fixed screening intervals and every women enters at 50
  my.Data <- data.frame(id = rep(1:N, each = Fup), N_crib = rep(1:Fup, N),
                        u = rep(u, each = Fup), Tentry = rep(Tenter, each = Fup),
                        Tsoj = rep(Tsoj, each = Fup), FP_1 = 0, FP = 0, cens = 0)
  rm(u, u_soj, Tsoj, Tenter)
  
  # dataset without left-truncation
  my.Data_50 <- modify_data(my.Data) %>%
    select(id, N_crib, time, time_cat, FP, state_un, state, select)
  rm(my.Data)
  
  ## MSM model
  msm_50[[.nsim]] <- model_msm_noFP(my.Data_50)
  
  if (.nsim %in% seq(0, Nsim - 1, 5)){
    save(msm_50, file = paste0('temp.',OFILE))
  }
  rm(my.Data_50)
  gc(reset = T)
}
save(msm_50, file = OFILE)