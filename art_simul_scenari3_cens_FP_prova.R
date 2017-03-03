# Simulation 4: Women with a time-depending FP as step function and random censoring
library(dplyr)
library(survival)
library(msm)

source("aux_fun.R")

FP_cens_creation <- function(.data, ncrib, data_ant, prob_cens_no, prob_cens_fp, prob_fp){
  adata <- .data %>% filter(N_crib == ncrib)
  adata$cens <- data_ant$cens # keep the previous censoring
  adata$FP_1 <- data_ant$FP_1 # keep the previous FP
  # Censoring for woman without FP
  adata$cens[adata$cens == 0 & adata$FP_1 == 0] <-
    with(adata, rbinom(sum(cens == 0 & FP_1 == 0), 1, prob_cens_no))
  # Censoring for woman with FP
  adata$cens[adata$cens == 0 & adata$FP_1 == 1] <-
    with(adata, rbinom(sum(cens == 0 & FP_1 == 1), 1, prob_cens_fp))
  # first FP
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
msm_50_pci <- list()
msm_50_covar <- list()
msm_left <- list() # list of MSM models with left-truncation
msm_left_pci <- list()
msm_left_covar <- list()
cox_50 <- list() # list of Cox models without left-truncation
cox_left <- list() # list of Cox models with left-truncation
disc_50 <- list() # list of Cox models without left-truncation
disc_left <- list() # list of Cox models with left-truncation
msm_50_4 <- list() # list of MSM models without left-truncation
msm_left_4 <- list() # list of MSM models with left-truncation
msm_left_4_pci <- list()
msm_left_4_covar <- list()
cox_50_4 <- list() # list of Cox models without left-truncation
cox_left_4 <- list() # list of Cox models with left-truncation
disc_50_4 <- list() # list of Cox models without left-truncation
disc_left_4 <- list() # list of Cox models with left-truncation
for (.nsim in 1:Nsim){
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
  # Creation of time-depending FP as step function
  # Cens and FP in the first observation
  my.Data_1 <- my.Data %>% filter(N_crib == 1)
  my.Data_1$cens <- 0 # No censoring in the first observation
  my.Data_1$FP_1 <- rbinom(N, 1, 9.610631/100) # probability from roman et al
  # Cens and FP in the second observation
  my.Data_2 <- FP_cens_creation(my.Data, 2, my.Data_1, 1 - 0.853, 1 - 0.793, 3.973549/100)
  # my.Data_2 <- my.Data %>% filter(N_crib == 2)
  # my.Data_2$cens <- my.Data_1$cens # keep the previous censoring
  # my.Data_2$FP_1 <- my.Data_1$FP_1 # keep the previous FP
  # # Censoring for woman without FP
  # my.Data_2$cens[my.Data_2$cens == 0 & my.Data_2$FP_1 == 0] <-
  #   rbinom(sum(my.Data_2$cens == 0 & my.Data_2$FP_1 == 0), 1, 1 - 0.853)
  # # Censoring for woman with FP
  # my.Data_2$cens[my.Data_2$cens == 0 & my.Data_2$FP_1 == 1] <-
  #   rbinom(sum(my.Data_2$cens == 0 & my.Data_2$FP_1 == 1), 1, 1 - 0.793)
  # first FP
  # my.Data_2$FP_1[my.Data_2$FP_1 == 0 & my.Data_2$cens == 0] <-
  #   rbinom(sum(my.Data_2$FP_1 == 0 & my.Data_2$cens == 0), 1, 3.973549/100)
  # Cens and FP for third observation
  my.Data_3 <- FP_cens_creation(my.Data, 3, my.Data_2, 1 - 0.888, 1 - 0.847, 3.057781/100)
  # my.Data_3 <- my.Data %>% filter(N_crib == 3)
  # my.Data_3$cens <- my.Data_2$cens
  # my.Data_3$FP_1 <- my.Data_2$FP_1
  # my.Data_3$cens[my.Data_3$cens == 0 & my.Data_3$FP_1 == 0] <-
  #   rbinom(sum(my.Data_3$cens == 0 & my.Data_3$FP_1 == 0), 1, 1 - 0.888)
  # my.Data_3$cens[my.Data_3$cens == 0 & my.Data_3$FP_1 == 1] <-
  #   rbinom(sum(my.Data_3$cens == 0 & my.Data_3$FP_1 == 1), 1, 1 - 0.847)
  # my.Data_3$FP_1[my.Data_3$FP_1 == 0 & my.Data_3$cens == 0] <-
  #   rbinom(sum(my.Data_3$FP_1 == 0 & my.Data_3$cens == 0), 1, 3.057781/100)
  # Cens and FP for forth observation
  my.Data_4 <- FP_cens_creation(my.Data, 4, my.Data_3, 1 - 0.899, 1 - 0.863, 2.798582/100)
  # my.Data_4 <- my.Data %>% filter(N_crib == 4)
  # my.Data_4$cens <- my.Data_3$cens
  # my.Data_4$FP_1 <- my.Data_3$FP_1
  # my.Data_4$cens[my.Data_4$cens == 0 & my.Data_4$FP_1 == 0] <-
  #   rbinom(sum(my.Data_4$cens == 0 & my.Data_4$FP_1 == 0), 1, 1 - 0.899)
  # my.Data_4$cens[my.Data_4$cens == 0 & my.Data_4$FP_1 == 1] <-
  #   rbinom(sum(my.Data_4$cens == 0 & my.Data_4$FP_1 == 1), 1, 1 - 0.863)
  # my.Data_4$FP_1[my.Data_4$FP_1 == 0 & my.Data_4$cens == 0] <-
  #   rbinom(sum(my.Data_4$FP_1 == 0 & my.Data_4$cens == 0), 1, 2.798582/100)
  my.Data_5 <- FP_cens_creation(my.Data, 5, my.Data_4, 1 - 0.920, 1 - 0.897, 2.410195/100)
  my.Data_6 <- FP_cens_creation(my.Data, 6, my.Data_5, 1 - 0.936, 1 - 0.918, 2.410195/100)
  my.Data_7 <- FP_cens_creation(my.Data, 7, my.Data_6, 1 - 0.960, 1 - 0.946, 2.410195/100)
  my.Data_8 <- FP_cens_creation(my.Data, 8, my.Data_7, 1 - 0.960, 1 - 0.946, 2.410195/100)
  my.Data_9 <- FP_cens_creation(my.Data, 9, my.Data_8, 1 - 0.960, 1 - 0.946, 2.410195/100)
  my.Data_10 <- FP_cens_creation(my.Data, 10, my.Data_9, 1 - 0.960, 1 - 0.946, 2.410195/100)
  # new dataset with censoring, first FP and lag of first FP
  my.Data <- rbind(my.Data_1, my.Data_2, my.Data_3, my.Data_4, my.Data_5, my.Data_6, my.Data_7,
                   my.Data_8, my.Data_9, my.Data_10) %>% arrange(id, N_crib) %>%
    mutate(FP = ifelse(N_crib == 1, 0, lag(FP_1))) %>% filter(cens == 0)
  
  # dataset without left-truncation
  my.Data_50 <- modify_data(my.Data)
  # dataset with left-truncation
  my.Data_left <- modify_data(my.Data, T)
  
  ## MSM model
  msm_50[[.nsim]] <- model_msm(my.Data_50, pci = F)
  msm_left[[.nsim]] <- model_msm(my.Data_left, pci = F)
  msm_50_4[[.nsim]] <- model_msm(my.Data_50 %>% filter(N_crib <= 4), pci = F)
  msm_left_4[[.nsim]] <- model_msm(my.Data_left %>% filter(N_crib <= 4), pci = F)
  msm_50_pci[[.nsim]] <- model_msm(my.Data_50)
  msm_left_pci[[.nsim]] <- model_msm(my.Data_left)
  msm_left_4_pci[[.nsim]] <- model_msm(my.Data_left %>% filter(N_crib <= 4))
  msm_50_covar[[.nsim]] <- model_msm(my.Data_50, covar = T)
  msm_left_covar[[.nsim]] <- model_msm(my.Data_left, covar = T)
  msm_left_4_covar[[.nsim]] <- model_msm(my.Data_left %>% filter(N_crib <= 4), covar = T)
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
save(msm_50, msm_50_pci, msm_50_covar, msm_left, msm_left_pci, msm_left_covar,
     msm_50_4, msm_left_4, msm_left_4_pci, msm_left_4_covar,
     cox_50, cox_left, cox_50_4, cox_left_4,
     disc_50, disc_left, disc_50_4, disc_left_4,
     file = OFILE)