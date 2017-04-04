#Afegir sensibilitat 85%
library(dplyr)
library(msm)

FP_creation <- function(.data, ncrib, data_ant, prob_fp){
  adata <- subset(.data, N_crib == ncrib)
  adata$FP_1 <- data_ant$FP_1 # keep the previous FP
  sel_FP <- adata$FP_1 == 0
  adata$FP_1[sel_FP] <- rbinom(sum(sel_FP), 1, prob_fp)
  adata
}

#Nsim <- 500 # number of simulation runs
Nsim <- 100 # number of simulation runs
N <- 62000 # number of women in the study
Fup <- 11
left <- F # sense left-truncation

# scale_param <- 688.1962
# shape_param <- 1
exp1_param <- (1/scale_param)^shape_param
exp2_param <- 1/4 #paràmetre exponencial
prob_fp_1 <- 9.610631/100
prob_fp_2 <- 3.973549/100
prob_fp_3 <- 3.057781/100
prob_fp_4 <- 2.798582/100
prob_fp_5 <- 2.410195/100

#msm_50 <- list() # list of MSM models without left-truncation

res_com <- list()
res_obs <- list()

# set.seed(n_seed[.nsim])
for (i in 1:Nsim){
  # Uniform for exponential distribution (Tpre)
  u <- runif(N)
  # sojourn time exponential distribution
  u_soj <- runif(N)
  Tsoj <- -log(u_soj)/exp2_param
  Tenter <- 0
  
  # create data set
  # assume fixed screening intervals and every women enters at 50
  my.Data <- data.frame(id = rep(1:N, each = Fup),
                        N_crib = rep(1:Fup, N),
                        Tentry = Tenter,
                        u = rep(u, each = Fup),
                        # Tpre = rep(Tpre, each = Fup),
                        Tsoj = rep(Tsoj, each = Fup),
                        obs_type = 1)
  rm(u, u_soj, Tsoj, Tenter)
  
  ## time at visits
  # my.Data$time <- with(my.Data, (N_crib - 1)*2) # time at visit
  my.Data$time <- with(my.Data, Tentry*left + (N_crib - 1)*2)

  if(fp_simul){
    # Creation of time-depending FP as step function
    my.Data_1 <- subset(my.Data, N_crib == 1)
    my.Data_1$FP_1 <- rbinom(N, 1, prob_fp_1) # probability from roman et al
    my.Data_2 <- FP_creation(my.Data, 2, my.Data_1, prob_fp_2)
    my.Data_3 <- FP_creation(my.Data, 3, my.Data_2, prob_fp_3)
    my.Data_4 <- FP_creation(my.Data, 4, my.Data_3, prob_fp_4)
    my.Data_5 <- FP_creation(my.Data, 5, my.Data_4, prob_fp_5)
    my.Data_6 <- FP_creation(my.Data, 6, my.Data_5, prob_fp_5)
    my.Data_7 <- FP_creation(my.Data, 7, my.Data_6, prob_fp_5)
    my.Data_8 <- FP_creation(my.Data, 8, my.Data_7, prob_fp_5)
    my.Data_9 <- FP_creation(my.Data, 9, my.Data_8, prob_fp_5)
    my.Data_10 <- FP_creation(my.Data, 10, my.Data_9, prob_fp_5)
    my.Data <- rbind(my.Data_1, my.Data_2, my.Data_3, my.Data_4, my.Data_5, my.Data_6, my.Data_7,
                     my.Data_8, my.Data_9, my.Data_10)
    rm(my.Data_1, my.Data_2, my.Data_3, my.Data_4, my.Data_5, my.Data_6, my.Data_7, my.Data_8,
       my.Data_9, my.Data_10)
    # my.Data <- my.Data[order(my.Data$id, my.Data$N_crib),]
    # my.Data$FP[2:dim(my.Data)[1]] <- my.Data$FP_1[1:(dim(my.Data)[1] - 1)]
    # my.Data$FP[my.Data$N_crib == 1] <- 0
    # # time of FP (No sé com fer-ho en base)
    # my.Data$time_FP <- with(my.Data, min(1000*(FP_1 == 0) + (FP_1 == 1)*time))
    # Ho reescric en dplyr syntaxi
    my.Data <- my.Data %>% arrange(id, N_crib)
    # Calculem el FP i el seu temps. Separem per cada id
    my.Data <- my.Data %>% group_by(id) %>%
      mutate(FP = ifelse(N_crib == 1, 0, lag(FP_1)),
             time_FP = min(1000*(FP_1 == 0) + (FP_1 == 1)*time))
  }
  
  # time of pre-clinical entrance
  if (trans12_exp){
    my.Data$Tpre <- -log(my.Data$u)/exp1_param
  } else{
    my.Data$Tpre <- (-log(my.Data$u)/exp1_param)^(1/shape_param)
    if(fp_simul){
      sel_FP <- with(my.Data, -log(u) >= exp1_param*time_FP^shape_param)
      my.Data$Tpre[sel_FP] <- with(my.Data[sel_FP,], ((-log(u) - exp1_param*time_FP^shape_param +
                                                         exp1_param*HR_FP*time_FP^shape_param)/
                                                        (exp1_param*HR_FP))^(1/shape_param))
    }
  }
  
  # time of clinical entrance
  my.Data$Tclin <- with(my.Data, Tpre + Tsoj)
  
  #Eliminem els prevalents
  my.Data <- subset(my.Data, !(Tpre <= Tentry))
  
  ## remove screening observations after Tclin+2
  ## hence Tclin < 6 included in data, even if not observed
  ## in order to see how estimates perform if follow-up continued after screen detected cancer
  my.Data <- subset(my.Data, time <= (Tclin+2))
  
  ## state variable at visit; 0, 1 or 2: no cancer, pre-clinical cancer, and clinical cancer
  my.Data$state_un <- with(my.Data, as.numeric(Tpre < time) + as.numeric(Tclin < time))
  my.Data$state <- my.Data$state_un
  
  id_fora <- c()
  my.Data$Tpre2 <- my.Data$Tpre
  if(sensibilitat){
    for (in_crib in 2:11){
      # Select the SCD for screening "in_crib" and without SCD
      sel_pre <- my.Data$N_crib == in_crib & my.Data$state_un == 1 & !(my.Data$id %in% id_fora)
      # Calculate the state (observed) with a sensibility of 0.85
      my.Data$state[sel_pre] <- rbinom(sum(sel_pre), 1, 0.85)
      sel_pre_id <- my.Data$id[sel_pre][my.Data$state[sel_pre] == 0]
      my.Data$Tpre2[my.Data$id %in% sel_pre_id] <-
        my.Data$Tpre2[my.Data$id %in% sel_pre_id] + 2
      id_fora <- c(id_fora, my.Data$id[my.Data$N_crib == in_crib & my.Data$state == 1])
      # new SCD database
    }
  }
  
  ## by choosing time > Tpre+2, the interval cancers have select=1
  # my.Data$observat <- with(my.Data, time <= (Tpre + 2))
  my.Data$observat <- with(my.Data, time <= (Tpre2 + 2))
  
  ## adjust time in last row if it refers to clinical cancer
  sel_clin <- my.Data$state == 2
  my.Data$time[sel_clin] <- my.Data$Tclin[sel_clin]
  my.Data$obs_type[sel_clin] <- 3
  
  #Cut time at 20 years
  my.Data <- subset(my.Data, time <= 20)
  
  # amy.Data <- my.Data
  
  if (afegir_CC){
    ## afegim la observació CC
    my.Data_CC <- subset(my.Data, N_crib == 1 & Tpre < 20)
    my.Data_CC$obs_type <- 3
    my.Data_CC$time <- my.Data_CC$Tpre
    my.Data_CC$state <- 1
    my.Data_CC$observat <- F
    my.Data <- rbind(my.Data, my.Data_CC)
    my.Data <- my.Data[order(my.Data$id, my.Data$time),]
  }
  
  # ## check whether data are ok
  # ## show first couple of individuals that left initial state
  # tmp <- subset(my.Data,state>0)[1:20,]
  # subset(my.Data, id %in% tmp$id)
  # ## table of observed states
  # table(my.Data$state)
  # ## table of number of visits per individual
  # table(table(my.Data$id))
  
  ###### include only de select==1 rows
  my.Data$state <- my.Data$state + 1
  my.Data$state_un <- my.Data$state_un + 1  
  
  ##### Initial values transition matrix
  Q <- rbind ( c(-0.004, 0.004, 0), 
               c(0,   -0.55, 0.55),
               c(0,    0,   0) )
  # Q.crude <- crudeinits.msm(state ~ time, subject = id, data = my.Data, qmatrix = Q)
  if (pci){
    if (fixed_pci){
      if(fp_simul){
        mod1.com <- msm(formula = state_un ~ time, subject = id, data = my.Data, qmatrix = Q,
                        obstype = obs_type, covariates = ~ FP, pci = c(4, 8, 12, 16),
                        fixedpars = c(6, 8, 10, 12),
                        method = "BFGS", control = list(fnscale = 10000, maxit = 10000))
        # sojourn.msm(mod1.com)
        mod1.obs <- msm(formula = state ~ time, subject = id, data = subset(my.Data, observat),
                        qmatrix = Q, obstype = obs_type, covariates = ~ FP,
                        pci = c(4, 8, 12, 16), fixedpars = c(6, 8, 10, 12), #Si fixem el FP 
                        method = "BFGS", control = list(fnscale = 10000, maxit = 10000))
        # sojourn.msm(mod1.com)
      } else{
        if(fp_simul){
          mod1.com <- msm(formula = state_un ~ time, subject = id, data = my.Data, qmatrix = Q,
                          obstype = obs_type, covariates = ~ FP, pci = c(4, 8, 12, 16),
                          fixedpars = c(4, 6, 8, 10),
                          method = "BFGS", control = list(fnscale = 10000, maxit = 10000))
          # sojourn.msm(mod1.com)
          mod1.obs <- msm(formula = state ~ time, subject = id, data = subset(my.Data, observat),
                          qmatrix = Q, obstype = obs_type, covariates = ~ FP,
                          pci = c(4, 8, 12, 16), fixedpars = c(4, 6, 8, 10),
                          method = "BFGS", control = list(fnscale = 10000, maxit = 10000))
          # sojourn.msm(mod1.com)
        } else{
          mod1.com <- msm(formula = state_un ~ time, subject = id, data = my.Data, qmatrix = Q,
                          obstype = obs_type, pci = c(4, 8, 12, 16), fixedpars = c(4, 6, 8, 10),
                          method = "BFGS", control = list(fnscale = 10000, maxit = 10000))
          # sojourn.msm(mod1.com)
          mod1.obs <- msm(formula = state ~ time, subject = id, data = subset(my.Data, observat),
                          qmatrix = Q, obstype = obs_type, pci = c(4, 8, 12, 16),
                          fixedpars = c(4, 6, 8, 10),
                          method = "BFGS", control = list(fnscale = 10000, maxit = 10000))
          # sojourn.msm(mod1.com)
        }
      }
    } else{
      if(fp_simul){
        mod1.com <- msm(formula = state_un ~ time, subject = id, data = my.Data, qmatrix = Q,
                        obstype = obs_type, covariates = ~ FP, pci = c(4, 8, 12, 16),
                        method = "BFGS", control = list(fnscale = 10000, maxit = 10000))
        # sojourn.msm(mod1.com)
        mod1.obs <- msm(formula = state ~ time, subject = id, data = subset(my.Data, observat),
                        qmatrix = Q, obstype = obs_type, covariates = ~ FP, pci = c(4, 8, 12, 16),
                        method = "BFGS", control = list(fnscale = 10000, maxit = 10000))
        # sojourn.msm(mod1.com)
      } else{
        mod1.com <- msm(formula = state_un ~ time, subject = id, data = my.Data, qmatrix = Q,
                        obstype = obs_type, pci = c(4, 8, 12, 16),
                        method = "BFGS", control = list(fnscale = 10000, maxit = 10000))
        # sojourn.msm(mod1.com)
        mod1.obs <- msm(formula = state ~ time, subject = id, data = subset(my.Data, observat),
                        qmatrix = Q, obstype = obs_type, pci = c(4, 8, 12, 16),
                        method = "BFGS", control = list(fnscale = 10000, maxit = 10000))
      }
    }
  } else{
    mod1.com <- msm(state_un ~ time, subject = id, data = my.Data, qmatrix = Q, obstype = obs_type,
                    method = "BFGS", control = list(fnscale = 10000, maxit = 10000))
    # sojourn.msm(mod1.com)
    mod1.obs <- msm(state ~ time, subject = id, data = subset(my.Data, observat), qmatrix = Q,
                    obstype = obs_type, method = "BFGS",
                    control = list(fnscale = 10000, maxit = 10000))
    # sojourn.msm(mod1.com)
  }
  
  res_com[[i]] <- mod1.com$estimates.t #results calculations
  res_obs[[i]] <- mod1.obs$estimates.t #results calculations
}

save(res_com, res_obs, file = OUT_FILE)
