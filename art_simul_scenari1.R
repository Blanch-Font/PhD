# Simulation 1: Women without FP
library(dplyr)
# Basic dplyr syntaxi:
#  - dplyr provides the %>% operator. x %>% f(y) turns into f(x, y).
#  - mutate: function to create new variables
#  - filter: function to remove rows that not satisfy a condition
#  - group_by: apply next functions to each categories of a variable
#  - ungroup: remove the function group_by
#  - select: keep/remove variables (columns)
#  - arrange: order the database according to variable
#  - summarise: It collapses a data frame to a single row
#  - bind_rows: similar to rbind
library(msm)
source("aux_fun.R")

# Calculating the scale parameter for the simulation model
aux_scale <- (1/scale_param)^shape_param

msm_50 <- list() # list of MSM models without left-truncation
msm_left <- list() # list of MSM models with left-truncation
msm_50_4 <- list() # list of MSM models without left-truncation and 4 obs
msm_left_4 <- list() # list of MSM models with left-truncation and 4 obs
for (.nsim in 1:Nsim){
  # time to preclinical cancer weibull distribution
  u <- runif(N)
  Tpre <- (-log(u)/aux_scale)^(1/shape_param)
  # sojourn time exponential distribution
  u_soj <- runif(N)
  Tsoj <- -log(u_soj)/0.25
  # Clinical time = Time to pre-clinical + sojourn time
  Tclin <- Tpre + Tsoj
  # for late entry, Tenter <- 0 if all enter at 50
  Tenter <- runif(N, 0, 15)
  # maximum number of visits
  Fup <- 10
  # create data set
  my.Data <- data.frame(id = rep(1:N, each = Fup), N_crib = rep(1:Fup, N),
                        Tentry = rep(Tenter, each = Fup), Tpre = rep(Tpre, each = Fup),
                        Tclin = rep(Tclin, each = Fup), Tsoj = rep(Tsoj, each = Fup))
  # dataset without left-truncation
  my.Data_50 <- modify_data_1(my.Data)
  # dataset with left-truncation
  my.Data_left <- modify_data_1(my.Data, T)
  
  ## MSM model
  msm_50[[.nsim]] <- model_msm_noFP(my.Data_50)
  msm_left[[.nsim]] <- model_msm_noFP(my.Data_left)
  msm_50_4[[.nsim]] <- model_msm_noFP(my.Data_50 %>% filter(N_crib <= 4), pci = F)
  msm_left_4[[.nsim]] <- model_msm_noFP(my.Data_left %>% filter(N_crib <= 4))
}
save(msm_50, msm_left, msm_50_4, msm_left_4, file = OFILE)