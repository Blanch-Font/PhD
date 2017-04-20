#Afegir sensibilitat 85%
library(dplyr)
library(msm)
library(survival)

FP_cens_creation <- function(.data, ncrib, data_ant, prob_cens, prob_fp){
  adata <- subset(.data, N_crib == ncrib)
  adata$cens <- data_ant$cens # keep the previous censoring
  adata$FP_1 <- data_ant$FP_1 # keep the previous FP
  if (length(prob_cens) == 1){
    adata$cens[adata$cens == 0] <- rbinom(N - sum(data_ant$cens), 1, prob_cens)
  } else {
    # Censoring for woman without FP
    adata$cens[adata$cens == 0 & adata$FP_1 == 0] <-
      with(adata, rbinom(sum(cens == 0 & FP_1 == 0), 1, prob_cens[1]))
    # Censoring for woman with FP
    adata$cens[adata$cens == 0 & adata$FP_1 == 1] <-
      with(adata, rbinom(sum(cens == 0 & FP_1 == 1), 1, prob_cens[2]))
  }
  # first FP
  adata$FP_1[adata$FP_1 == 0 & adata$cens == 0] <-
    with(adata, rbinom(sum(FP_1 == 0 & cens == 0), 1, prob_fp))
  adata
}

# Discret-time Model
model_discret <- function(.data){
  # Multinomial
  .data_un <- .data %>% filter(state_un != 2) %>% group_by(id, N_crib) %>%
    summarise(type = 'scd', FP = max(FP), event = max(as.numeric(state_un == 1))) %>%
    bind_rows(.data %>% group_by(id, N_crib) %>%
                summarise(type = 'ic', FP = max(FP), event = max(as.numeric(state_un == 2))))
  .data_sel <- .data %>% filter(observat & state != 2) %>% group_by(id, N_crib) %>%
    summarise(type = 'scd', FP = max(FP), event = max(as.numeric(state == 1))) %>%
    bind_rows(.data %>% filter(observat) %>% group_by(id, N_crib) %>%
                summarise(type = 'ic', FP = max(FP), event = max(as.numeric(state == 2))))
  mod_al <- glm(formula = event ~ type:(as.factor(N_crib) + FP - 1),
                family = binomial(link = "logit"), data = .data_un)
  mod_ob <- glm(formula = event ~ type:(as.factor(N_crib) + FP - 1),
                family = binomial(link = "logit"), data = .data_sel)
  #mod_al_cloglog <- glm(formula = event ~ type:(as.factor(N_crib) + FP - 1),
  #                      family = binomial(link = "cloglog"), data = .data_un)
  #mod_ob_cloglog <- glm(formula = event ~ type:(as.factor(N_crib) + FP - 1),
  #                      family = binomial(link = "cloglog"), data = .data_sel)
  
  # 2-state model
  # model without extending the observation of IC
  .data_no <- .data %>% filter(observat) %>%
    mutate(N_crib = ifelse(state == 2, as.integer(lag(N_crib)), N_crib))
  mod_no <- glm(formula = I(state > 0) ~ as.factor(N_crib) + FP - 1,
                family = binomial(link = "logit"), data = .data_no)
  #mod_no_cloglog <- glm(formula = I(state > 0) ~ as.factor(N_crib) + FP - 1,
  #                      family = binomial(link = "cloglog"), data = .data_no)
  # model with extending the observation of IC
  .data_si <- .data %>% filter(observat) %>%
    mutate(N_crib = ifelse(state == 2, as.integer(lag(N_crib) + 1), N_crib))
  mod_si <- glm(formula = I(state > 0) ~ as.factor(N_crib) + FP - 1,
                family = binomial(link = "logit"), data = .data_si)
  #mod_si_cloglog <- glm(formula = I(state > 0) ~ as.factor(N_crib) + FP - 1,
  #                      family = binomial(link = "cloglog"), data = .data_si)
  list(mod_al = list(coef = coef(mod_al), confint = confint(mod_al)),
       mod_ob = list(coef = coef(mod_ob), confint = confint(mod_ob)),
       mod_no = list(coef = coef(mod_no), confint = confint(mod_no)),
       mod_si = list(coef = coef(mod_si), confint = confint(mod_si)))
  #mod_al_cloglog = list(coef = coef(mod_al_cloglog), confint = confint(mod_al_cloglog)),
  #mod_ob_cloglog = list(coef = coef(mod_ob_cloglog), confint = confint(mod_ob_cloglog)),
  #mod_no_cloglog = list(coef = coef(mod_no_cloglog), confint = confint(mod_no_cloglog)),
  #mod_si_cloglog = list(coef = coef(mod_si_cloglog), confint = confint(mod_si_cloglog)))
  # list(mod_al = list(coef = coef(mod_al)),
  #      mod_ob = list(coef = coef(mod_ob)),
  #      mod_no = list(coef = coef(mod_no)),
  #      mod_si = list(coef = coef(mod_si)))
}

#Nsim <- 500 # number of simulation runs
# Nsim <- 100 # number of simulation runs
# N <- 62000 # number of women in the study
Fup <- 11
# left <- F # sense left-truncation

# scale_param <- 688.1962
# shape_param <- 1
exp1_param <- (1/scale_param)^shape_param
exp2_param <- 1/4 #paràmetre exponencial

msm_mod <- list() # list of MSM models without left-truncation
msm_mod_4 <- list() # list of MSM models without left-truncation
cox_mod <- list()
cox_mod_4 <- list()
disc_mod <- list()
disc_mod_4 <- list()

# set.seed(n_seed[.nsim])
for (.nsim in 1:Nsim){
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
    my.Data_1$cens <- 0 # No censoring in the first observation
    my.Data_1$FP_1 <- rbinom(N, 1, prob_fp_1) # probability from roman et al
    my.Data_2 <- FP_cens_creation(my.Data, 2, my.Data_1, prob_cens2, prob_fp_2)
    my.Data_3 <- FP_cens_creation(my.Data, 3, my.Data_2, prob_cens3, prob_fp_3)
    my.Data_4 <- FP_cens_creation(my.Data, 4, my.Data_3, prob_cens4, prob_fp_4)
    my.Data_5 <- FP_cens_creation(my.Data, 5, my.Data_4, prob_cens5, prob_fp_5)
    my.Data_6 <- FP_cens_creation(my.Data, 6, my.Data_5, prob_cens6, prob_fp_5)
    my.Data_7 <- FP_cens_creation(my.Data, 7, my.Data_6, prob_cens7, prob_fp_5)
    my.Data_8 <- FP_cens_creation(my.Data, 8, my.Data_7, prob_cens7, prob_fp_5)
    my.Data_9 <- FP_cens_creation(my.Data, 9, my.Data_8, prob_cens7, prob_fp_5)
    my.Data_10 <- FP_cens_creation(my.Data, 10, my.Data_9, prob_cens7, prob_fp_5)
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
  
  ## by choosing time > Tpre+2, the interval cancers have observat=1
  # my.Data$observat <- with(my.Data, time <= (Tpre + 2))
  my.Data$observat <- with(my.Data, time <= (Tpre2 + 2))
  
  ## adjust time in last row if it refers to clinical cancer
  sel_clin <- my.Data$state == 2
  my.Data$time[sel_clin] <- my.Data$Tclin[sel_clin]
  my.Data$obs_type[sel_clin] <- 3
  
  #Cut time at 20 years
  my.Data <- subset(my.Data, time <= 20)

  #For IC, redefine N_crib as N_crib - 1
  sel_IC <- my.Data$state_un == 2
  my.Data$N_crib[sel_IC] <- my.Data$N_crib[sel_IC] - 1
  
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
  
  ## Discret models
  disc_mod[[.nsim]] <- model_discret(my.Data)
  disc_mod_4[[.nsim]] <- model_discret(my.Data %>% filter(N_crib <= 4))

  if (.nsim %in% seq(0, Nsim - 1, 5)){
    save(disc_mod, disc_mod_4, file = paste0('temp.',OFILE))
  }
  rm(my.Data)
  gc(reset = T)
}

save(disc_mod, disc_mod_4, file = OFILE)
