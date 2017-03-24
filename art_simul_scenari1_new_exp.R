# Simulation 4: Women with a time-depending FP as step function and random censoring
library(dplyr)
library(survival)
library(msm)

source("aux_fun.R")

# Function to modify database
modify_data_new <- function(.data, left = F){
  # time at visit
  .data$time <- with(.data, Tentry*left + (N_crib - 1)*2)
  # time to first FP (if no FP -> time of FP = 1000)
  .data$time_FP <- with(.data, min(1000*(FP_1 == 0) + (FP_1 == 1)*time))
  # time of preclinical entrance (as Weibull) and no FP
  .data$Tpre <- with(.data, (-log(u)/aux_scale)^(1/shape_param))
  # modify the preclinical entrance depending on the time of FP
  sel <- with(.data, -log(u) < aux_scale*(time_FP)^shape_param)
  .data$Tpre[sel] <- with(.data[sel, ], ((-log(u) - aux_scale*(time_FP)^shape_param +
                                            aux_scale*HR_FP*(time_FP)^shape_param)/
                                           (aux_scale*HR_FP))^(1/shape_param))
  # time of clinical entrance
  .data$Tclin <- with(.data, Tpre + Tsoj)
  # state at visit (unobserved)
  .data$state_un <- with(.data, as.numeric(time > Tpre) + as.numeric(time > Tclin))
  # Separate the database in the three states (No, SCD and IC)
  # Database for No cancer
  .data_No <- .data %>% filter(state_un == 0)
  # state (observed) = No cancer
  .data_No$state <- 0
  # Database for SCD
  .data_SCD <- .data[.data$state_un == 1, ]
  .data_SCD_gen <- .data_SCD
  .data_SCD_gen$state <- 1
  # keep the first observation with a pre-clinical stage
  # id_fora <- c()
  # .data_SCD_gen <- NULL
  # for (i in 1:10){
  #   # Select the SCD for screening "i" and without SCD
  #   .adata_SCD <- .data_SCD[.data_SCD$N_crib == i & !(.data_SCD$id %in% id_fora), ]
  #   # Calculate the state (observed) with a sensibility of 0.85
  #   .adata_SCD$state <- rbinom(dim(.adata_SCD)[1], 1, 0.85)
  #   id_fora <- c(id_fora, .adata_SCD$id[.adata_SCD$state == 1])
  #   # new SCD database
  #   .data_SCD_gen <- .data_SCD_gen %>% bind_rows(.adata_SCD)
  # }
  # Database for IC
  .data_IC <- .data %>% filter(state_un > 0) %>%
    semi_join(.data %>% filter(state_un > 0) %>% group_by(id, state_un) %>%
                summarise(N_crib = min(N_crib)), by = c('id', 'N_crib')) %>%
    group_by(id) %>%
    summarise(N_crib = max(N_crib), u = first(u), Tentry = first(Tentry), Tsoj = first(Tsoj),
              cens = first(cens), FP_1 = first(FP_1), FP = first(FP), time = first(Tclin),
              time_FP = first(time_FP), Tpre = first(Tpre), Tclin = first(Tclin),
              state_un = 2, state = 2)
  # join in a one dataset
  .data <- bind_rows(.data_No, bind_rows(.data_SCD_gen, .data_IC)) %>% arrange(id, time)
  # keep screenings less than 70 years old (20 years in the study)
  .data <- .data %>% filter(time < 20)
  # Remove the prevalent cancer
  sel_can <- with(.data, (state == lag(state)) & (state > 0))
  sel_can[.data$N_crib == 1] <- F
  .data <- .data[!sel_can, ]
  # 
  sel_sel <- with(.data, (Tpre + 2 < time) | (state_un == 2 & lag(state_un) == 1))
  .data$select <- T
  .data$select[sel_sel] <- F
  if (left) .data <- .data[.data$Tentry < .data$Tpre,] #eliminem
  # # Add one observation for Preclinical stage
  # sel_preclin <- with(.data, (0 < state_un) & (lag(state_un) == 0))
  # sel_preclin[.data$N_crib == 1] <- F
  # .data_scd <- .data[sel_preclin, ]
  # .data_scd$time <- .data_scd$Tpre
  # .data_scd$state_un <- 1
  # .data_scd$state <- 0
  # .data_scd$select <- F
  # # Return a dataset with preclinical state and order by id and time
  # .data <- bind_rows(.data, .data_scd) %>% arrange(id, time)
  .data$time_cat <- cut(x = .data$time, breaks = c(0, 5, 10, 15, 25), right = FALSE)
  .data
}

# Parameters
aux_scale <- (1/scale_param)^shape_param
Fup <- 10 # maximum number of visits
HR_FP <- 1

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
  my.Data_50 <- modify_data_new(my.Data) %>%
    select(id, N_crib, time, time_cat, FP, state_un, state, select)
  rm(my.Data)
  
  ## MSM model
  msm_50[[.nsim]] <- model_msm_noFP(my.Data_50, pci = F)
  
  if (.nsim %in% seq(0, Nsim - 1, 5)){
    save(msm_50, file = paste0('temp.',OFILE))
  }
  rm(my.Data_50)
  gc(reset = T)
}
save(msm_50, file = OFILE)
