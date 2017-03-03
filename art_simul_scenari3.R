# Simulation 3: Women with a time-depending FP as step function
library(dplyr)
library(survival)
library(msm)

source("aux_fun.R")
FP_creation <- function(.data, ncrib, data_ant, prob){
  adata <- .data %>% filter(N_crib == ncrib)
  adata$FP_1 <- data_ant$FP_1 # keep the previuos FP
  # for free-FP woman, calculate the first FP
  adata$FP_1[adata$FP_1 == 0] <- rbinom(N - sum(data_ant$FP_1), 1, prob)
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
  # FP in the first observation
  my.Data_1 <- my.Data %>% filter(N_crib == 1)
  my.Data_1$FP_1 <- rbinom(N, 1, 9.610631/100) # probability from our data
  # FP in the second observation
  my.Data_2 <- FP_creation(my.Data, 2, my.Data_1, 3.973549/100)
  # my.Data_2 <- my.Data %>% filter(N_crib == 2)
  # my.Data_2$FP_1 <- my.Data_1$FP_1 # keep the previuos FP
  # # for free-FP woman, calculate the first FP
  # my.Data_2$FP_1[my.Data_2$FP_1 == 0] <- rbinom(N - sum(my.Data_1$FP_1), 1, 3.973549/100)
  # FP in the third observation
  my.Data_3 <- FP_creation(my.Data, 3, my.Data_2, 3.057781/100)
  # my.Data_3 <- my.Data %>% filter(N_crib == 3)
  # my.Data_3$FP_1 <- my.Data_2$FP_1
  # my.Data_3$FP_1[my.Data_3$FP_1 == 0] <- rbinom(N - sum(my.Data_2$FP_1), 1, 3.057781/100)
  # FP in the forth observation
  my.Data_4 <- FP_creation(my.Data, 4, my.Data_3, 2.798582/100)
  # my.Data_4 <- my.Data %>% filter(N_crib == 4)
  # my.Data_4$FP_1 <- my.Data_3$FP_1
  # my.Data_4$FP_1[my.Data_4$FP_1 == 0] <- rbinom(N - sum(my.Data_4$FP_1), 1, 2.798582/100)
  # FALTA Generar de 5 a 10 observacions
  my.Data_5 <- FP_creation(my.Data, 5, my.Data_4, 2.410195/100)
  my.Data_6 <- FP_creation(my.Data, 6, my.Data_5, 2.410195/100)
  my.Data_7 <- FP_creation(my.Data, 7, my.Data_6, 2.410195/100)
  my.Data_8 <- FP_creation(my.Data, 8, my.Data_7, 2.410195/100)
  my.Data_9 <- FP_creation(my.Data, 9, my.Data_8, 2.410195/100)
  my.Data_10 <- FP_creation(my.Data, 10, my.Data_9, 2.410195/100)
  
  # new dataset with first FP and lag of first FP
  my.Data <- rbind(my.Data_1, my.Data_2, my.Data_3, my.Data_4, my.Data_5, my.Data_6, my.Data_7,
                   my.Data_8, my.Data_9, my.Data_10) %>% arrange(id, N_crib) %>%
    mutate(FP = ifelse(N_crib == 1, 0, lag(FP_1)))
  
  # dataset without left-truncation
  my.Data_50 <- modify_data(my.Data)
  # dataset with left-truncation
  my.Data_left <- modify_data(my.Data, T)
  
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
     disc_50, disc_left, disc_50_4, disc_left_4, file = OFILE)