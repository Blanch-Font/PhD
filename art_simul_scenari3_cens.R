# Simulation 4: Women with a time-depending FP as step function and random censoring
library(dplyr)
library(survival)
library(msm)

source("aux_fun.R")

FP_cens_creation <- function(.data, ncrib, data_ant, prob_cens, prob_fp){
  adata <- .data %>% filter(N_crib == ncrib)
  adata$cens <- data_ant$cens # keep censoring from last observation
  adata$cens[adata$cens == 0] <- rbinom(N - sum(data_ant$cens), 1, prob_cens) 
  adata$FP_1 <- data_ant$FP_1 # keep the previous FP
  adata$FP_1[adata$FP_1 == 0 & adata$cens == 0] <- 
    with(adata, rbinom(sum(FP_1 == 0 & cens == 0), 1, prob_fp))
  adata
}

# Parameters
# Nsim <- 1000 # number of simulation runs
# N <- 15000 # number of women in the study
# scale_param <- 0.002 # scale parameter for Weibull distribution
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
                        Tsoj = rep(Tsoj, each = Fup))
  rm(u, u_soj, Tsoj, Tenter)
  # Creation of time-depending FP as step function
  # Cens and FP in the first observation
  my.Data_1 <- my.Data %>% filter(N_crib == 1)
  my.Data_1$cens <- 0 # No censoring in the first observation
  my.Data_1$FP_1 <- rbinom(N, 1, 9.610631/100) # probability from roman et al
  # Cens and FP in the second observation
  my.Data_2 <- FP_cens_creation(my.Data, 2, my.Data_1, 1 - 0.817, 3.973549/100)
  # my.Data_2 <- my.Data %>% filter(N_crib == 2)
  # my.Data_2$cens <- my.Data_1$cens # keep censoring from last observation
  # my.Data_2$cens[my.Data_2$cens == 0] <- rbinom(N - sum(my.Data_1$cens), 1, 1 - 0.817) 
  # my.Data_2$FP_1 <- my.Data_1$FP_1 # keep the previous FP
  # my.Data_2$FP_1[my.Data_2$FP_1 == 0 & my.Data_2$cens == 0] <-
  #   rbinom(sum(my.Data_2$FP_1 == 0 & my.Data_2$cens == 0), 1, 3.973549/100)
  # Cens and FP in the third observation
  my.Data_3 <- FP_cens_creation(my.Data, 3, my.Data_2, 1 - 0.866, 3.057781/100)
  # my.Data_3 <- my.Data %>% filter(N_crib == 3)
  # my.Data_3$cens <- my.Data_2$cens
  # my.Data_3$cens[my.Data_3$cens == 0] <- rbinom(N - sum(my.Data_2$cens), 1, 1 - 0.866)
  # my.Data_3$FP_1 <- my.Data_2$FP_1
  # my.Data_3$FP_1[my.Data_3$FP_1 == 0 & my.Data_3$cens == 0] <-
  #   rbinom(sum(my.Data_3$FP_1 == 0 & my.Data_3$cens == 0), 1, 3.057781/100)
  # Cens and FP in the forth obseration
  my.Data_4 <- FP_cens_creation(my.Data, 4, my.Data_3, 1 - 0.881, 2.798582/100)
  # my.Data_4 <- my.Data %>% filter(N_crib == 4)
  # my.Data_4$cens <- my.Data_3$cens
  # my.Data_4$cens[my.Data_4$cens == 0] <- rbinom(N - sum(my.Data_3$cens), 1, 1 - 0.881)
  # my.Data_4$FP_1 <- my.Data_3$FP_1
  # my.Data_4$FP_1[my.Data_4$FP_1 == 0 & my.Data_4$cens == 0] <-
  #   rbinom(sum(my.Data_4$FP_1 == 0 & my.Data_4$cens == 0), 1, 2.798582/100)
  my.Data_5 <- FP_cens_creation(my.Data, 5, my.Data_4, 1 - 0.905, 2.410195/100)
  my.Data_6 <- FP_cens_creation(my.Data, 6, my.Data_5, 1 - 0.925, 2.410195/100)
  my.Data_7 <- FP_cens_creation(my.Data, 7, my.Data_6, 1 - 0.956, 2.410195/100)
  my.Data_8 <- FP_cens_creation(my.Data, 8, my.Data_7, 1 - 0.956, 2.410195/100)
  my.Data_9 <- FP_cens_creation(my.Data, 9, my.Data_8, 1 - 0.956, 2.410195/100)
  my.Data_10 <- FP_cens_creation(my.Data, 10, my.Data_9, 1 - 0.956, 2.410195/100)
  # new dataset with censoring, first FP and lag of first FP
  my.Data <- rbind(my.Data_1, my.Data_2, my.Data_3, my.Data_4, my.Data_5, my.Data_6, my.Data_7,
                   my.Data_8, my.Data_9, my.Data_10) %>% arrange(id, N_crib) %>%
    mutate(FP = ifelse(N_crib == 1, 0, lag(FP_1))) %>% filter(cens == 0)
  rm(my.Data_1, my.Data_2, my.Data_3, my.Data_4, my.Data_5, my.Data_6, my.Data_7, my.Data_8,
     my.Data_9, my.Data_10)
  
  # dataset without left-truncation
  my.Data_50 <- modify_data(my.Data) %>%
    select(id, N_crib, time, time_cat, FP, state_un, state, select)
  # dataset with left-truncation
  my.Data_left <- modify_data(my.Data, T) %>%
    select(id, N_crib, time, time_cat, FP, state_un, state, select)
  rm(my.Data)
  
  ## MSM model
  msm_50[[.nsim]] <- model_msm_ini(my.Data_50)
  msm_left[[.nsim]] <- model_msm_ini(my.Data_left)
  msm_50_4[[.nsim]] <- model_msm_ini(my.Data_50 %>% filter(N_crib <= 4), pci = F)
  msm_left_4[[.nsim]] <- model_msm_ini(my.Data_left %>% filter(N_crib <= 4))
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
  
  if (.nsim %in% seq(0, Nsim, 5)){
    save(msm_50, msm_left, cox_50, cox_left,
         msm_50_4, msm_left_4, cox_50_4, cox_left_4,
         disc_50, disc_left, disc_50_4, disc_left_4,
         file = paste0('temp.',OFILE))
  }
  rm(my.Data_50, my.Data_left)
  gc(reset = T)
}
save(msm_50, msm_left, cox_50, cox_left,
     msm_50_4, msm_left_4, cox_50_4, cox_left_4,
     disc_50, disc_left, disc_50_4, disc_left_4, file = OFILE)