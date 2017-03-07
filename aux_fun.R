# Function to modify database
modify_data_1 <- function(.data, left = F){
  # Calculate time and state at visit
  .data <- .data %>% 
    mutate(time = Tentry*left + (N_crib - 1)*2, # time at visit
           state_un = as.numeric(time > Tpre) + as.numeric(time > Tclin))
  # state = as.numeric(time > Tpre) + as.numeric(time > Tclin)) # state at visit
  # # remove screening observations after Tclin+2
  # # hence Tclin < 6 included in data, even if not observed
  # # in order to see how estimates perform if follow-up continued after SCD
  # .data <- .data %>% filter(time <= (Tclin + 2))
  .data_No <- .data %>% filter(state_un == 0)
  .data_No$state <- 0
  .data_SCD <- .data %>% filter(state_un == 1)
  id_fora <- c()
  .data_SCD_gen <- NULL
  for (i in 1:10){
    .adata_SCD <- .data_SCD %>% filter(N_crib == i & !(id %in% id_fora))
    .adata_SCD$state <- rbinom(dim(.adata_SCD)[1], 1, 0.85)
    id_fora <- c(id_fora, .adata_SCD$id[.adata_SCD$state == 1])
    .data_SCD_gen <- .data_SCD_gen %>% bind_rows(.adata_SCD)
  }
  #ens quedem amb el 1r registre amb IC
  .data_IC <- .data %>% filter(state_un > 0) %>%
    semi_join(.data %>% filter(state_un > 0) %>% group_by(id, state_un) %>%
                summarise(N_crib = min(N_crib)), by = c('id', 'N_crib')) %>% group_by(id) %>%
    summarise(N_crib = max(N_crib), Tentry = first(Tentry), Tpre = first(Tpre),
              Tclin = first(Tclin), Tsoj = first(Tsoj), time = Tclin, state_un = 2, state = 2)
  # .data_IC <- .data %>% filter(state_un == 2) %>% 
  #   semi_join(.data %>% filter(state_un == 2) %>% group_by(id) %>%
  #               summarise(N_crib = min(N_crib)), by = c('id', 'N_crib')) 
  # .data_IC$state <- 2
  .data <- .data_No %>% bind_rows(.data_SCD_gen) %>% bind_rows(.data_IC) %>%
    arrange(id, time) %>% group_by(id)
  # keep screenings less than 70 years old (20 years in the study)
  .data <- .data %>% filter(time < 20)
  # Remove the prevalent cancer
  .data <- .data %>%# group_by(id) %>%
    mutate(sel_can = ifelse(N_crib == 1, F, (state == lag(state)) & (state > 0))) %>%
    filter(!sel_can) %>% select(-sel_can) %>% ungroup()
  # column to see which rows are not observed (select=0)
  # by choosing time > Tpre+2, the interval cancers have select=1
  # .data$select <- T
  # .data[(.data$Tpre + 2) < .data$time, "select"] <- F
  .data <- .data %>%
    mutate(select = ifelse(Tpre + 2 < time | (state_un == 2 & lag(state_un) == 1), F, T))
  if (left) .data <- .data %>% filter(Tentry < Tpre) #eliminem
  # ## adjust time in last row if it refers to clinical cancer
  # .data[.data$Tclin < .data$time, "time"] <- .data[.data$Tclin < .data$time, "Tclin"]
  # Add one observation for Preclinical stage
  .data_scd <- .data %>% group_by(id) %>%
    mutate(sel_preclin = ifelse(N_crib == 1, F, (0 < state_un) & (lag(state_un) == 0))) %>%
    filter(sel_preclin) %>% select(-sel_preclin) %>%
    mutate(time = Tpre, state_un = 1, state = 0, select = F)
  # Return a dataset with preclinical state and order by id and time
  .data <- .data %>% bind_rows(.data_scd) %>% group_by(id) %>% arrange(id, time)
  .data$time_cat <- cut(x = .data$time, breaks = c(0, 5, 10, 15, 25), right = FALSE)
  .data
}

modify_data_2 <- function(.data, left = F){
  # Calculate time and state at visit
  .data <- .data %>% 
    mutate(time = Tentry*left + (N_crib - 1)*2, # time at visit
           state_un = as.numeric(time > Tpre) + as.numeric(time > Tclin))
  # state = as.numeric(time > Tpre) + as.numeric(time > Tclin)) # state at visit
  # # remove screening observations after Tclin+2
  # # hence Tclin < 6 included in data, even if not observed
  # # in order to see how estimates perform if follow-up continued after SCD
  # .data <- .data %>% filter(time <= (Tclin + 2))
  .data_No <- .data %>% filter(state_un == 0)
  .data_No$state <- 0
  .data_SCD <- .data %>% filter(state_un == 1)
  id_fora <- c()
  .data_SCD_gen <- NULL
  for (i in 1:10){
    .adata_SCD <- .data_SCD %>% filter(N_crib == i & !(id %in% id_fora))
    .adata_SCD$state <- rbinom(dim(.adata_SCD)[1], 1, 0.85)
    id_fora <- c(id_fora, .adata_SCD$id[.adata_SCD$state == 1])
    .data_SCD_gen <- .data_SCD_gen %>% bind_rows(.adata_SCD)
  }
  #ens quedem amb el 1r registre amb IC
  .data_IC <- .data %>% filter(state_un > 0) %>%
    semi_join(.data %>% filter(state_un > 0) %>% group_by(id, state_un) %>%
                summarise(N_crib = min(N_crib)), by = c('id', 'N_crib')) %>% group_by(id) %>%
    summarise(N_crib = max(N_crib), FP = first(FP), Tentry = first(Tentry), Tpre = first(Tpre),
              Tclin = first(Tclin), Tsoj = first(Tsoj), time = Tclin, state_un = 2, state = 2)
  # .data_IC <- .data %>% filter(state_un == 2) %>% 
  #   semi_join(.data %>% filter(state_un == 2) %>% group_by(id) %>%
  #               summarise(N_crib = min(N_crib)), by = c('id', 'N_crib')) 
  # .data_IC$state <- 2
  .data <- .data_No %>% bind_rows(.data_SCD_gen) %>% bind_rows(.data_IC) %>%
    arrange(id, time) %>% group_by(id)
  # keep screenings less than 70 years old (20 years in the study)
  .data <- .data %>% filter(time < 20)
  # Remove the prevalent cancer
  .data <- .data %>%# group_by(id) %>%
    mutate(sel_can = ifelse(N_crib == 1, F, (state == lag(state)) & (state > 0))) %>%
    filter(!sel_can) %>% select(-sel_can) %>% ungroup()
  # column to see which rows are not observed (select=0)
  # by choosing time > Tpre+2, the interval cancers have select=1
  # .data$select <- T
  # .data[(.data$Tpre + 2) < .data$time, "select"] <- F
  .data <- .data %>%
    mutate(select = ifelse(Tpre + 2 < time | (state_un == 2 & lag(state_un) == 1), F, T))
  if (left) .data <- .data %>% filter(Tentry < Tpre) #eliminem
  # ## adjust time in last row if it refers to clinical cancer
  # .data[.data$Tclin < .data$time, "time"] <- .data[.data$Tclin < .data$time, "Tclin"]
  # Add one observation for Preclinical stage
  .data_scd <- .data %>% group_by(id) %>%
    mutate(sel_preclin = ifelse(N_crib == 1, F, (0 < state_un) & (lag(state_un) == 0))) %>%
    filter(sel_preclin) %>% select(-sel_preclin) %>%
    mutate(time = Tpre, state_un = 1, state = 0, select = F)
  # Return a dataset with preclinical state and order by id and time
  .data <- .data %>% bind_rows(.data_scd) %>% group_by(id) %>% arrange(id, time)
  .data$time_cat <- cut(x = .data$time, breaks = c(0, 5, 10, 15, 25), right = FALSE)
  .data
}

# Function to modify database
modify_data <- function(.data, left = F){
  # Calculate time at visit, time of the first FP for each woman, 
  # time to pre-clinical as Austin (Stat in Med. DOI: 10.1002/sim.5452),
  # Clinical time = Time to pre-clinical + sojourn time
  #.data$time <- with(.data, Tentry*left + (N_crib - 1)*2) # time at visit
  #.data$time_FP <- with(.data, min(1000*(FP_1 == 0) + (FP_1 == 1)*time)) # time to first FP
  #.data$Tpre <- with(.data, ((-log(u) - aux_scale*time_FP^shape_param +
  #                           aux_scale*HR_FP*time_FP^shape_param)/
  #                          (aux_scale*HR_FP))^(1/shape_param)))
  #.data$Tpre[-log(.data$u) < aux_scale*.data$time_FP^shape_param] <-
  #  with(.data, (-log(u)/aux_scale)^(1/shape_param))
  #.data$Tclin <- with(.data, Tpre + Tsoj)
  .data <- .data %>% group_by(id) %>%
    mutate(time = Tentry*left + (N_crib - 1)*2, # time at visit
           time_FP = min(1000*(FP_1 == 0) + (FP_1 == 1)*time), # time to first FP
           Tpre = ifelse(-log(u)<aux_scale*time_FP^shape_param,
                         (-log(u)/aux_scale)^(1/shape_param),
                         ((-log(u) - aux_scale*time_FP^shape_param +
                             aux_scale*HR_FP*time_FP^shape_param)/
                            (aux_scale*HR_FP))^(1/shape_param)),
           Tclin = Tpre + Tsoj)
  # state at visit
  .data$state_un <- with(.data, as.numeric(time > Tpre) + as.numeric(time > Tclin))
  # # state according the past states
  # .data <- .data %>% group_by(id) %>% mutate(state = cummax(state))
  # keep all the observations without cancer
  # data_s0 <- .data %>% filter(state == 0)
  .data_No <- .data %>% filter(state_un == 0)
  .data_No$state <- 0
  .data_SCD <- .data %>% filter(state_un == 1)
  # keep the first observation with a pre-clinical stage
  # data_s1 <- .data %>% filter(state == 1) %>% group_by(id, state) %>%
  #   filter(time <= (min(Tpre) + 2))
  id_fora <- c()
  .data_SCD_gen <- NULL
  for (i in 1:10){
    .adata_SCD <- .data_SCD %>% filter(N_crib == i & !(id %in% id_fora))
    .adata_SCD$state <- rbinom(dim(.adata_SCD)[1], 1, 0.85)
    id_fora <- c(id_fora, .adata_SCD$id[.adata_SCD$state == 1])
    .data_SCD_gen <- .data_SCD_gen %>% bind_rows(.adata_SCD)
  }
  # keep the first observation with a clinical stage
  # data_s2 <- .data %>% filter(state == 2) %>% group_by(id, state) %>%
  #   filter(time <= (min(Tclin) + 2))
  .data_IC <- .data %>% filter(state_un > 0) %>%
    semi_join(.data %>% filter(state_un > 0) %>% group_by(id, state_un) %>%
                summarise(N_crib = min(N_crib)), by = c('id', 'N_crib')) %>% group_by(id) %>%
    summarise(N_crib = max(N_crib), u = first(u), Tentry = first(Tentry), Tsoj = first(Tsoj),
              cens = first(cens), FP_1 = first(FP_1), FP = first(FP), time = first(Tclin),
              time_FP = first(time_FP), Tpre = first(Tpre), Tclin = first(Tclin),
              state_un = 2, state = 2)
  # join in a one dataset
  # .data <- data_s0 %>% bind_rows(data_s1) %>% bind_rows(data_s2) %>% arrange(id, time)
  .data <- .data_No %>% bind_rows(.data_SCD_gen) %>% bind_rows(.data_IC) %>%
    arrange(id, time) %>% group_by(id)
  # keep screenings less than 70 years old (20 years in the study)
  .data <- .data %>% filter(time < 20)
  # Remove the prevalent cancer
  .data <- .data %>%# group_by(id) %>%
    mutate(sel_can = ifelse(N_crib == 1, F, (state == lag(state)) & (state > 0))) %>%
    filter(!sel_can) %>% select(-sel_can) %>% ungroup()
  # column to see which rows are not observed (select=0)
  # by choosing time > Tpre+2, the interval cancers have select=1
  # .data$select <- T
  # .data[(.data$Tpre + 2) < .data$time, "select"] <- F
  .data <- .data %>%
    mutate(select = ifelse(Tpre + 2 < time | (state_un == 2 & lag(state_un) == 1), F, T))
  if (left) .data <- .data %>% filter(Tentry < Tpre) #eliminem
  # # For left-truncation, remove the prevalent cancers
  # if (left) .data[.data$Tpre < .data$Tentry, "select"] <- F
  # # adjust time in last row if it refers to clinical cancer
  # .data[.data$Tclin < .data$time, "time"] <- .data[.data$Tclin < .data$time, "Tclin"]
  # # Add one observation for Preclinical stage
  .data_scd <- .data %>% group_by(id) %>%
    mutate(sel_preclin = ifelse(N_crib == 1, F, (0 < state_un) & (lag(state_un) == 0))) %>%
    filter(sel_preclin) %>% select(-sel_preclin) %>%
    mutate(time = Tpre, state_un = 1, state = 0, select = F)
  # Return a dataset with preclinical state and order by id and time
  .data <- .data %>% bind_rows(.data_scd) %>% group_by(id) %>% arrange(id, time)
  .data$time_cat <- cut(x = .data$time, breaks = c(0, 5, 10, 15, 25), right = FALSE)
  .data
}

# MSM models
model_msm_noFP <- function(.data, pci = T){
  # 3-state models
  # Modify state to satisfy MSM assumptions
  .data_msm <- .data %>% group_by(id) %>% mutate(event_un = state_un + 1, event = state + 1)
  qmat3 <- rbind(c(0.9, 0.1, 0), c(0, 0.6, 0.4), c(0, 0, 0))
  # 2-state model
  # model without extending the observation of IC
  .data_no <- .data %>% filter(select) %>% mutate(event = as.numeric(state > 0) + 1)
  qmat2 <- rbind(c(0.9, 0.1), c(0, 0))
  # model with extending the observation of IC
  # new time for IC corresponding the last visit + 2 yr
  .data_si <- .data %>% filter(select) %>%
    mutate(time = ifelse(state == 2, lag(time) + 2, time), event = as.numeric(state > 0) + 1)
  if (pci){
    mod_al <- msm(formula = event_un ~ time, subject = id, data = .data_msm, qmatrix = qmat3,
                  deathexact = 3, # covinits = list(FP = c(log(HR_FP), 0)),
                  pci = c(4, 8, 12, 16), fixedpars = c(4, 6, 8, 10),
                  control = list(fnscale = 100, reltol = 1e-7))
    mod_ob <- msm(formula = event ~ time, subject = id, data = .data_msm %>% filter(select),
                  qmatrix = qmat3, deathexact = 3, pci = c(4, 8, 12, 16),
                  fixedpars = c(4, 6, 8, 10), control = list(fnscale = 100, reltol = 1e-7))
    mod_no <- msm(formula = event ~ time, subject = id, data = .data_no, qmatrix = qmat2,
                  pci = c(4, 8, 12, 16), deathexact = 2)
    mod_si <- msm(formula = event ~ time, subject = id, data = .data_si, qmatrix = qmat2,
                  pci = c(4, 8, 12, 16), deathexact = 2)
  } else{
    mod_al <- msm(formula = event_un ~ time, subject = id, data = .data_msm, qmatrix = qmat3,
                  deathexact = 3, control = list(fnscale = 100, reltol = 1e-7))
    mod_ob <- msm(formula = event ~ time, subject = id, data = .data_msm %>% filter(select),
                  qmatrix = qmat3, deathexact = 3, control = list(fnscale = 100, reltol = 1e-7))
    mod_no <- msm(formula = event ~ time, subject = id, data = .data_no, qmatrix = qmat2,
                  deathexact = 2)
    mod_si <- msm(formula = event ~ time, subject = id, data = .data_si, qmatrix = qmat2,
                  deathexact = 2)
  }
  # Results:
  #  - HZ: hazard ratios
  #  - qtrans: transition intensity matrix for each period
  list(mod_al = list(HZ = hazard.msm(mod_al),
                     qtrans = list(t1 = qmatrix.msm(mod_al,
                                                    covariates = list(timeperiod = "[-Inf,4)")),
                                   t2 = qmatrix.msm(mod_al,
                                                    covariates = list(timeperiod = "[4,8)")),
                                   t3 = qmatrix.msm(mod_al,
                                                    covariates = list(timeperiod = "[8,12)")),
                                   t4 = qmatrix.msm(mod_al,
                                                    covariates = list(timeperiod = "[12,16)")),
                                   t5 = qmatrix.msm(mod_al,
                                                    covariates = list(timeperiod = "[16,Inf)")))),
       mod_ob = list(HZ = hazard.msm(mod_ob),
                     qtrans = list(t1 = qmatrix.msm(mod_ob,
                                                    covariates = list(timeperiod = "[-Inf,4)")),
                                   t2 = qmatrix.msm(mod_ob,
                                                    covariates = list(timeperiod = "[4,8)")),
                                   t3 = qmatrix.msm(mod_ob,
                                                    covariates = list(timeperiod = "[8,12)")),
                                   t4 = qmatrix.msm(mod_ob,
                                                    covariates = list(timeperiod = "[12,16)")),
                                   t5 = qmatrix.msm(mod_ob,
                                                    covariates = list(timeperiod = "[16,Inf)")))),
       mod_no = list(HZ = hazard.msm(mod_no),
                     qtrans = list(t1 = qmatrix.msm(mod_no,
                                                    covariates = list(timeperiod = "[-Inf,4)")),
                                   t2 = qmatrix.msm(mod_no,
                                                    covariates = list(timeperiod = "[4,8)")),
                                   t3 = qmatrix.msm(mod_no,
                                                    covariates = list(timeperiod = "[8,12)")),
                                   t4 = qmatrix.msm(mod_no,
                                                    covariates = list(timeperiod = "[12,16)")),
                                   t5 = qmatrix.msm(mod_no,
                                                    covariates = list(timeperiod = "[16,Inf)")))),
       mod_si = list(HZ = hazard.msm(mod_si),
                     qtrans = list(t1 = qmatrix.msm(mod_si,
                                                    covariates = list(timeperiod = "[-Inf,4)")),
                                   t2 = qmatrix.msm(mod_si,
                                                    covariates = list(timeperiod = "[4,8)")),
                                   t3 = qmatrix.msm(mod_si,
                                                    covariates = list(timeperiod = "[8,12)")),
                                   t4 = qmatrix.msm(mod_si,
                                                    covariates = list(timeperiod = "[12,16)")),
                                   t5 = qmatrix.msm(mod_si,
                                                    covariates = list(timeperiod = "[16,Inf)")))))
}
# MSM model
model_msm <- function(.data, pci = T, covar = F){
  # 3-state models
  # Modify state to satisfy MSM assumptions
  .data_msm <- .data %>% group_by(id) %>% mutate(event_un = state_un + 1, event = state + 1)
  qmat3 <- rbind(c(0.9, 0.1, 0), c(0, 0.6, 0.4), c(0, 0, 0))
  # 2-state model
  # model without extending the observation of IC
  .data_no <- .data %>% filter(select) %>% mutate(event = as.numeric(state > 0) + 1)
  qmat2 <- rbind(c(0.9, 0.1), c(0, 0))
  # model with extending the observation of IC
  # new time for IC corresponding the last visit + 2 yr
  .data_si <- .data %>% filter(select) %>%
    mutate(time = ifelse(state == 2, lag(time) + 2, time), event = as.numeric(state > 0) + 1)
  
  # with all data
  if (pci){
    if (covar){
      mod_al <- msm(formula = event_un ~ time, subject = id, data = .data_msm, qmatrix = qmat3,
                    deathexact = 3, covariates = ~ FP + time_cat,
                    control = list(fnscale = 100, reltol = 1e-7))
      mod_ob <- msm(formula = event ~ time, subject = id, data = .data_msm %>% filter(select),
                    qmatrix = qmat3, deathexact = 3, covariates = ~ FP + time_cat,
                    control = list(fnscale = 100, reltol = 1e-7))
      mod_no <- msm(formula = event ~ time, subject = id, data = .data_no, qmatrix = qmat2,
                    covariates = ~ FP + time_cat, deathexact = 2)
      mod_si <- msm(formula = event ~ time, subject = id, data = .data_si, qmatrix = qmat2,
                    covariates = ~ FP + time_cat, deathexact = 2)
    } else{
      mod_al <- msm(formula = event_un ~ time, subject = id, data = .data_msm, qmatrix = qmat3,
                    deathexact = 3, covariates = ~ FP, pci = c(5, 10, 15),
                    control = list(fnscale = 100, reltol = 1e-7))
      mod_ob <- msm(formula = event ~ time, subject = id, data = .data_msm %>% filter(select),
                    qmatrix = qmat3, deathexact = 3, covariates = ~ FP, pci = c(5, 10, 15),
                    control = list(fnscale = 100, reltol = 1e-7))
      mod_no <- msm(formula = event ~ time, subject = id, data = .data_no, qmatrix = qmat2,
                    covariates = ~ FP, pci = c(5, 10, 15), deathexact = 2)
      mod_si <- msm(formula = event ~ time, subject = id, data = .data_si, qmatrix = qmat2,
                    covariates = ~ FP, pci = c(5, 10, 15), deathexact = 2)
    }
  } else{
    mod_al <- msm(formula = event_un ~ time, subject = id, data = .data_msm, qmatrix = qmat3,
                  deathexact = 3, covariates = ~ FP, 
                  control = list(fnscale = 100, reltol = 1e-7))
    mod_ob <- msm(formula = event ~ time, subject = id, data = .data_msm %>% filter(select),
                  qmatrix = qmat3, deathexact = 3, covariates = ~ FP,
                  control = list(fnscale = 100, reltol = 1e-7))
    mod_no <- msm(formula = event ~ time, subject = id, data = .data_no, qmatrix = qmat2,
                  covariates = ~ FP, deathexact = 2)
    mod_si <- msm(formula = event ~ time, subject = id, data = .data_si, qmatrix = qmat2,
                  covariates = ~ FP, deathexact = 2)
  }

  list(mod_al = list(HZ = hazard.msm(mod_al),
                     qtrans = list(t1 = qmatrix.msm(mod_al,
                                                    covariates = list(timeperiod = "[-Inf,5)")),
                                   t2 = qmatrix.msm(mod_al,
                                                    covariates = list(timeperiod = "[5,10)")),
                                   t3 = qmatrix.msm(mod_al,
                                                    covariates = list(timeperiod = "[10,15)")),
                                   t4 = qmatrix.msm(mod_al,
                                                    covariates = list(timeperiod = "[15,Inf)")))),
       mod_ob = list(HZ = hazard.msm(mod_ob),
                     qtrans = list(t1 = qmatrix.msm(mod_ob,
                                                    covariates = list(timeperiod = "[-Inf,5)")),
                                   t2 = qmatrix.msm(mod_ob,
                                                    covariates = list(timeperiod = "[5,10)")),
                                   t3 = qmatrix.msm(mod_ob,
                                                    covariates = list(timeperiod = "[10,15)")),
                                   t4 = qmatrix.msm(mod_ob,
                                                    covariates = list(timeperiod = "[15,Inf)")))),
       mod_no = list(HZ = hazard.msm(mod_no),
                     qtrans = list(t1 = qmatrix.msm(mod_no,
                                                    covariates = list(timeperiod = "[-Inf,5)")),
                                   t2 = qmatrix.msm(mod_no,
                                                    covariates = list(timeperiod = "[5,10)")),
                                   t3 = qmatrix.msm(mod_no,
                                                    covariates = list(timeperiod = "[10,15)")),
                                   t4 = qmatrix.msm(mod_no,
                                                    covariates = list(timeperiod = "[15,Inf)")))),
       mod_si = list(HZ = hazard.msm(mod_si),
                     qtrans = list(t1 = qmatrix.msm(mod_si,
                                                    covariates = list(timeperiod = "[-Inf,5)")),
                                   t2 = qmatrix.msm(mod_si,
                                                    covariates = list(timeperiod = "[5,10)")),
                                   t3 = qmatrix.msm(mod_si,
                                                    covariates = list(timeperiod = "[10,15)")),
                                   t4 = qmatrix.msm(mod_si,
                                                    covariates = list(timeperiod = "[15,Inf)")))))
}

# MSM model amb valors inicials
model_msm_ini <- function(.data, pci = TRUE){
  # 3-state models
  # Modify state to satisfy MSM assumptions
  .data_msm <- .data %>% group_by(id) %>% mutate(event_un = state_un + 1, event = state + 1)
  qmat3 <- rbind(c(0.9, 0.1, 0), c(0, 0.6, 0.4), c(0, 0, 0))
  # 2-state model
  # model without extending the observation of IC
  .data_no <- .data %>% filter(select) %>% mutate(event = as.numeric(state > 0) + 1)
  qmat2 <- rbind(c(0.9, 0.1), c(0, 0))
  # model with extending the observation of IC
  # new time for IC corresponding the last visit + 2 yr
  .data_si <- .data %>% filter(select) %>%
    mutate(time = ifelse(state == 2, lag(time) + 2, time), event = as.numeric(state > 0) + 1)
  
  # with all data
  if (pci){
    mod_al <- msm(formula = event_un ~ time, subject = id, data = .data_msm, qmatrix = qmat3,
                  deathexact = 3, covariates = ~ FP,# covinits = list(FP = c(log(HR_FP), 0)),
                  pci = c(4, 8, 12, 16), fixedpars = c(6, 8, 10, 12),
                  control = list(fnscale = 100, reltol = 1e-7))
    mod_ob <- msm(formula = event ~ time, subject = id, data = .data_msm %>% filter(select),
                  qmatrix = qmat3, deathexact = 3, covariates = ~ FP,
                  pci = c(4, 8, 12, 16),# covinits = list(FP = c(log(HR_FP), 0)),
                  fixedpars = c(6, 8, 10, 12), control = list(fnscale = 100, reltol = 1e-7))
    mod_no <- msm(formula = event ~ time, subject = id, data = .data_no, qmatrix = qmat2,
                  covariates = ~ FP, #covinits = list(FP = log(HR_FP)),
                  pci = c(4, 8, 12, 16), deathexact = 2)
    mod_si <- msm(formula = event ~ time, subject = id, data = .data_si, qmatrix = qmat2,
                  covariates = ~ FP,# covinits = list(FP = log(HR_FP)),
                  pci = c(4, 8, 12, 16), deathexact = 2)
  } else{
    mod_al <- msm(formula = event_un ~ time, subject = id, data = .data_msm, qmatrix = qmat3,
                  deathexact = 3, covariates = ~ FP, 
                  control = list(fnscale = 100, reltol = 1e-7))
    mod_ob <- msm(formula = event ~ time, subject = id, data = .data_msm %>% filter(select),
                  qmatrix = qmat3, deathexact = 3, covariates = ~ FP,
                  control = list(fnscale = 100, reltol = 1e-7))
    mod_no <- msm(formula = event ~ time, subject = id, data = .data_no, qmatrix = qmat2,
                  covariates = ~ FP, deathexact = 2)
    mod_si <- msm(formula = event ~ time, subject = id, data = .data_si, qmatrix = qmat2,
                  covariates = ~ FP, deathexact = 2)
  }
  
  list(mod_al = list(HZ = hazard.msm(mod_al),
                     qtrans = list(t1 = qmatrix.msm(mod_al,
                                                    covariates = list(timeperiod = "[-Inf,4)")),
                                   t2 = qmatrix.msm(mod_al,
                                                    covariates = list(timeperiod = "[4,8)")),
                                   t3 = qmatrix.msm(mod_al,
                                                    covariates = list(timeperiod = "[8,12)")),
                                   t4 = qmatrix.msm(mod_al,
                                                    covariates = list(timeperiod = "[12,16)")),
                                   t5 = qmatrix.msm(mod_al,
                                                    covariates = list(timeperiod = "[16,Inf)")))),
       mod_ob = list(HZ = hazard.msm(mod_ob),
                     qtrans = list(t1 = qmatrix.msm(mod_ob,
                                                    covariates = list(timeperiod = "[-Inf,4)")),
                                   t2 = qmatrix.msm(mod_ob,
                                                    covariates = list(timeperiod = "[4,8)")),
                                   t3 = qmatrix.msm(mod_ob,
                                                    covariates = list(timeperiod = "[8,12)")),
                                   t4 = qmatrix.msm(mod_ob,
                                                    covariates = list(timeperiod = "[12,16)")),
                                   t5 = qmatrix.msm(mod_ob,
                                                    covariates = list(timeperiod = "[16,Inf)")))),
       mod_no = list(HZ = hazard.msm(mod_no),
                     qtrans = list(t1 = qmatrix.msm(mod_no,
                                                    covariates = list(timeperiod = "[-Inf,4)")),
                                   t2 = qmatrix.msm(mod_no,
                                                    covariates = list(timeperiod = "[4,8)")),
                                   t3 = qmatrix.msm(mod_no,
                                                    covariates = list(timeperiod = "[8,12)")),
                                   t4 = qmatrix.msm(mod_no,
                                                    covariates = list(timeperiod = "[12,16)")),
                                   t5 = qmatrix.msm(mod_no,
                                                    covariates = list(timeperiod = "[16,Inf)")))),
       mod_si = list(HZ = hazard.msm(mod_si),
                     qtrans = list(t1 = qmatrix.msm(mod_si,
                                                    covariates = list(timeperiod = "[-Inf,4)")),
                                   t2 = qmatrix.msm(mod_si,
                                                    covariates = list(timeperiod = "[4,8)")),
                                   t3 = qmatrix.msm(mod_si,
                                                    covariates = list(timeperiod = "[8,12)")),
                                   t4 = qmatrix.msm(mod_si,
                                                    covariates = list(timeperiod = "[12,16)")),
                                   t5 = qmatrix.msm(mod_si,
                                                    covariates = list(timeperiod = "[16,Inf)")))))
}

# MSM model amb valors inicials
model_msm_inca <- function(.data, pci = TRUE){
  # 3-state models
  # Modify state to satisfy MSM assumptions
  # .data_msm <- .data %>% mutate(event = as.numeric(tipuscan))
  .data_msm <- .data %>% mutate(event = as.numeric(tipuscan == "CC") + 1) %>%
    bind_rows(.data %>% filter(tipuscan == "CI") %>% mutate(T0 = T1, event = 3)) %>%
    arrange(Mujer.id, N.Crib3, T0)
  qmat3 <- rbind(c(0.9, 0.1, 0), c(0, 0.6, 0.4), c(0, 0, 0))
  # 2-state model
  # model without extending the observation of IC
  .data_no <- .data_msm %>% mutate(event = as.numeric(event > 1) + 1)
  qmat2 <- rbind(c(0.9, 0.1), c(0, 0))
  # model with extending the observation of IC
  # new time for IC corresponding the last visit + 2 yr
  .data_si <- .data_no %>% 
    mutate(T0 = ifelse(tipuscan == "CI" & event == 2, lag(T0) + 2, T0))
  
  # with all data
  if (pci){
    mod_3 <- msm(formula = event ~ T0, subject = Mujer.id, data = .data_msm, qmatrix = qmat3,
                 deathexact = 3, covariates = ~ FP.msm, pci = c(4, 8, 12, 16),
                 fixedpars = c(4, 6, 8, 10, 12), control = list(fnscale = 100, reltol = 1e-7))
    mod_no <- msm(formula = event ~ T0, subject = Mujer.id, data = .data_no, qmatrix = qmat2,
                  covariates = ~ FP.msm, pci = c(4, 8, 12, 16), deathexact = 2)
    mod_si <- msm(formula = event ~ T0, subject = Mujer.id, data = .data_si, qmatrix = qmat2,
                  covariates = ~ FP.msm, pci = c(4, 8, 12, 16), deathexact = 2)
  } else{
    mod_3 <- msm(formula = event ~ T0, subject = Mujer.id, data = .data_msm, qmatrix = qmat3,
                 deathexact = 3, covariates = ~ FP.msm,
                 control = list(fnscale = 100, reltol = 1e-7))
    mod_no <- msm(formula = event ~ T0, subject = id, data = .data_no, qmatrix = qmat2,
                  covariates = ~ FP.msm, deathexact = 2)
    mod_si <- msm(formula = event ~ T0, subject = Mujer.id, data = .data_si, qmatrix = qmat2,
                  covariates = ~ FP.msm, deathexact = 2)
  }
  
  list(mod_3 = list(HZ = hazard.msm(mod_3),
                    qtrans = list(t1 = qmatrix.msm(mod_3,
                                                   covariates = list(timeperiod = "[-Inf,4)")),
                                  t2 = qmatrix.msm(mod_3,
                                                   covariates = list(timeperiod = "[4,8)")),
                                  t3 = qmatrix.msm(mod_3,
                                                   covariates = list(timeperiod = "[8,12)")),
                                  t4 = qmatrix.msm(mod_3,
                                                   covariates = list(timeperiod = "[12,16)")),
                                  t5 = qmatrix.msm(mod_3,
                                                   covariates = list(timeperiod = "[16,Inf)")))),
       mod_no = list(HZ = hazard.msm(mod_no),
                     qtrans = list(t1 = qmatrix.msm(mod_no,
                                                    covariates = list(timeperiod = "[-Inf,4)")),
                                   t2 = qmatrix.msm(mod_no,
                                                    covariates = list(timeperiod = "[4,8)")),
                                   t3 = qmatrix.msm(mod_no,
                                                    covariates = list(timeperiod = "[8,12)")),
                                   t4 = qmatrix.msm(mod_no,
                                                    covariates = list(timeperiod = "[12,16)")),
                                   t5 = qmatrix.msm(mod_no,
                                                    covariates = list(timeperiod = "[16,Inf)")))),
       mod_si = list(HZ = hazard.msm(mod_si),
                     qtrans = list(t1 = qmatrix.msm(mod_si,
                                                    covariates = list(timeperiod = "[-Inf,4)")),
                                   t2 = qmatrix.msm(mod_si,
                                                    covariates = list(timeperiod = "[4,8)")),
                                   t3 = qmatrix.msm(mod_si,
                                                    covariates = list(timeperiod = "[8,12)")),
                                   t4 = qmatrix.msm(mod_si,
                                                    covariates = list(timeperiod = "[12,16)")),
                                   t5 = qmatrix.msm(mod_si,
                                                    covariates = list(timeperiod = "[16,Inf)")))))
}

# # Cox model
# model_cox <- function(.data){
#   # 3-state model
#   # t0: time to begin the observation
#   # t1: time to end the observation
#   # event: state at the end of observation
#   .data_cox <- .data %>% mutate(t0 = time, t1 = lead(time), event = lead(state)) %>% 
#     filter(!is.na(t1)) %>% mutate(estrat = ifelse(event == 2, "clin", "pre"))
#   # with all information
#   mod_al <- coxph(Surv(t0, t1, event > 0) ~ FP*strata(estrat == "pre"), data = .data_cox)
#   # hiding some information
#   mod_ob <- coxph(Surv(t0, t1, event > 0) ~ FP*strata(estrat == "pre"), data = .data_cox, 
#                   subset = select)
#   # 2-state model
#   # model without extending the observation of IC
#   .data_no <- .data %>% filter(select) %>%
#     mutate(t0 = time, t1 = lead(time), event = lead(state > 0)) %>% filter(!is.na(t1))
#   mod_no <- coxph(Surv(t0, t1, event, type = "counting") ~ FP, data = .data_no)
#   # model with extending the observation of IC
#   .data_si <- .data %>% filter(select) %>%
#     mutate(t0 = time, t1 = ifelse(lead(state) == 2, time + 2, lead(time)),
#            event = lead(state > 0)) %>% filter(!is.na(t1))
#   mod_si <- coxph(Surv(t0, t1, event, type = "counting") ~ FP, data = .data_si)
#   list(mod_al = list(coef = coef(mod_al), confint = confint(mod_al)), 
#        mod_ob = list(coef = coef(mod_ob), confint = confint(mod_ob)),
#        mod_no = list(coef = coef(mod_no), confint = confint(mod_no)),
#        mod_si = list(coef = coef(mod_si), confint = confint(mod_si)))
# }
# Cox Model
model_cox <- function(.data){
  # causa-específic
  .data_un <- .data %>% filter(state_un != 2) %>% group_by(id) %>%
    mutate(type = 'scd', t0 = time, t1 = lead(time), event = lead(state_un)) %>%
    bind_rows(.data %>% group_by(id) %>%
                mutate(type = 'ic', t0 = time, t1 = lead(time), event = lead(state_un) == 2)) %>%
    filter(!is.na(t1))
  .data_un <- .data_un %>% filter(event == 0) %>%
    bind_rows(.data_un %>% filter(event == 1) %>% group_by(type, id) %>%
                summarise(N_crib = first(N_crib),  time = min(time), time_cat = first(time_cat),
                          FP = first(FP), state_un = first(state_un), state = first(state),
                          select = first(select), t0 = min(t0), t1 = min(t1),
                          event = min(event))) %>%
    arrange(type, id, t0)
  .data_sel <- .data %>% filter(select & state != 2) %>% group_by(id) %>%
    mutate(type = 'scd', t0 = time, t1 = lead(time), event = lead(state)) %>%
    bind_rows(.data %>% filter(select) %>% group_by(id) %>%
                mutate(type = 'ic', t0 = time, t1 = lead(time), event = lead(state) == 2)) %>%
    filter(!is.na(t1))
  mod_al <- coxph(Surv(t0, t1, event == 1) ~ strata(type):FP + cluster(id), data = .data_un)
  mod_ob <- coxph(Surv(t0, t1, event == 1) ~ strata(type):FP + cluster(id), data = .data_sel)
  
  # 2-state model
  # model without extending the observation of IC
  .data_no <- .data %>% filter(select) %>% group_by(id) %>%
    mutate(t0 = time, t1 = lead(time), event = lead(state > 0)) %>% filter(!is.na(t1))
  mod_no <- coxph(Surv(t0, t1, event) ~ FP + cluster(id), data = .data_no)
  # model with extending the observation of IC
  .data_si <- .data %>% filter(select) %>% group_by(id) %>%
    mutate(t0 = time, t1 = ifelse(lead(state) == 2, time + 2, lead(time)),
           event = lead(state > 0)) %>% filter(!is.na(t1))
  mod_si <- coxph(Surv(t0, t1, event) ~ FP + cluster(id), data = .data_si)
  list(mod_al = list(coef = coef(mod_al), confint = confint(mod_al)),
       mod_ob = list(coef = coef(mod_ob), confint = confint(mod_ob)),
       mod_no = list(coef = coef(mod_no), confint = confint(mod_no)),
       mod_si = list(coef = coef(mod_si), confint = confint(mod_si)))
}

model_cox_inca <- function(.data){
  # causa-específic
  .data_cc <- .data %>% mutate(type = 'scd', event = tipuscan == "CC") %>%
    bind_rows(.data %>% mutate(type = 'ic', event = tipuscan == "CI"))
  mod_cc <- coxph(Surv(T0, T1, event) ~ strata(type):FP.msm + cluster(Mujer.id),
                  data = .data_cc)
  
  # 2-state model
  # model without extending the observation of IC
  .data_no <- .data
  .data_no$event <- .data_no$tipuscan != 'No'
  mod_no <- coxph(Surv(T0, T1, event) ~ FP.msm + cluster(Mujer.id), data = .data_no)
  # model with extending the observation of IC
  .data_si <- .data_no
  .data_si$T1[.data_si$tipuscan == 'CI'] <- .data_si$T0[.data_si$tipuscan == 'CI'] + 2
  mod_si <- coxph(Surv(T0, T1, event) ~ FP.msm + cluster(Mujer.id), data = .data_si)
  list(mod_cc = list(coef = coef(mod_cc), confint = confint(mod_cc)),
       mod_no = list(coef = coef(mod_no), confint = confint(mod_no)),
       mod_si = list(coef = coef(mod_si), confint = confint(mod_si)))
}

# Discret-time Model
model_discret <- function(.data){
  # Multinomial
  .data_un <- .data %>% filter(state_un != 2) %>% group_by(id, N_crib) %>%
    summarise(type = 'scd', FP = max(FP), event = max(as.numeric(state_un == 1))) %>%
    bind_rows(.data %>% group_by(id, N_crib) %>%
                summarise(type = 'ic', FP = max(FP), event = max(as.numeric(state_un == 2))))
  .data_sel <- .data %>% filter(select & state != 2) %>% group_by(id, N_crib) %>%
    summarise(type = 'scd', FP = max(FP), event = max(as.numeric(state == 1))) %>%
    bind_rows(.data %>% filter(select) %>% group_by(id, N_crib) %>%
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
  .data_no <- .data %>% filter(select) %>%
    mutate(N_crib = ifelse(state == 2, as.integer(lag(N_crib)), N_crib))
  mod_no <- glm(formula = I(state > 0) ~ as.factor(N_crib) + FP - 1,
                family = binomial(link = "logit"), data = .data_no)
  #mod_no_cloglog <- glm(formula = I(state > 0) ~ as.factor(N_crib) + FP - 1,
  #                      family = binomial(link = "cloglog"), data = .data_no)
  # model with extending the observation of IC
  .data_si <- .data %>% filter(select) %>%
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

model_discret_inca <- function(.data){
  # Multinomial
  .data_cc <- .data %>% mutate(type = 'scd', event = as.numeric(tipuscan == "CC")) %>%
    bind_rows(.data %>% mutate(type = 'ic', event = as.numeric(tipuscan == "CI")))
  .data_cc$N.Crib3[.data_cc$N.Crib3 > 4] <- 4
  mod_cc <- glm(formula = event ~ type:(as.factor(N.Crib3) + FP.msm - 1),
                family = binomial(link = "logit"), data = .data_cc)
  
  # 2-state model
  # model without extending the observation of IC
  .data_no <- .data
  .data_no$event <- .data_no$tipuscan != 'No'
  .data_si <- .data_no
  .data_si$event[.data_si$tipuscan == 'CI'] <- F
  .data_si <- .data_si %>%
    bind_rows(.data_si %>% filter(tipuscan == 'CI') %>%
                mutate(N.Crib3 = N.Crib3 + 1,  event = T)) %>% arrange(Mujer.id, N.Crib3)
  .data_no$N.Crib3[.data_no$N.Crib3 > 4] <- 4
  mod_no <- glm(formula = event ~ as.factor(N.Crib3) + FP.msm - 1,
                family = binomial(link = "logit"), data = .data_no)
  # model with extending the observation of IC
  mod_si <- glm(formula = event ~ as.factor(N.Crib3) + FP.msm - 1,
                family = binomial(link = "logit"), data = .data_si)
  list(mod_cc = list(coef = coef(mod_cc), confint = confint(mod_cc)),
       mod_no = list(coef = coef(mod_no), confint = confint(mod_no)),
       mod_si = list(coef = coef(mod_si), confint = confint(mod_si)))
}

print_qtrans <- function(models, pci, from, to, param){
  aux <- unlist(lapply(models, function(.mod) .mod$qtrans[[pci]]$estimates[from,to]))
  se <- sqrt(sum((aux - mean(aux))^2)/(length(aux) - 1))
  lo <- unlist(lapply(models, function(.mod) .mod$qtrans[[pci]]$L[from,to]))
  hi <- unlist(lapply(models, function(.mod) .mod$qtrans[[pci]]$U[from,to]))
  resul <- data.frame(mean = mean(aux), se, median = median(aux), 
                      p25 = quantile(aux, probs = 0.025), p975 = quantile(aux, probs = 0.975),
                      bias = mean(aux) - param, bias_tp = (mean(aux) - param)/param*100,
                      accuracy = (mean(aux) - param)^2 + se^2,
                      coverage = sum(lo < param & param <= hi)/length(lo)*100)
  row.names(resul) <- ""
  resul
}

print_qtrans_covar <- function(models, covar, from, to, param){
  if (covar == '[0,5)'){
    aux <- unlist(lapply(models, function(.mod) .mod$qtrans$t1$estimates[from,to]))
    lo <- unlist(lapply(models, function(.mod) .mod$qtrans$t1$L[from,to]))
    hi <- unlist(lapply(models, function(.mod) .mod$qtrans$t1$U[from,to]))
  } else{
    aux <- unlist(lapply(models, function(.mod){
      .mod$qtrans$t1$estimates[from,to]*.mod$HZ[[paste0('time_cat', covar)]]
    }))
    lo <- unlist(lapply(models, function(.mod){
      .mod$qtrans$t1$L[from,to]*.mod$HZ[[paste0('time_cat', covar)]]
    }))
    hi <- unlist(lapply(models, function(.mod){
      .mod$qtrans$t1$U[from,to]*.mod$HZ[[paste0('time_cat', covar)]]
    }))
  }
  se <- sqrt(sum((aux - mean(aux))^2)/(length(aux) - 1))
  
  resul <- data.frame(mean = mean(aux), se, median = median(aux), 
                      p25 = quantile(aux, probs = 0.025), p975 = quantile(aux, probs = 0.975),
                      bias = mean(aux) - param, bias_tp = (mean(aux) - param)/param*100,
                      accuracy = (mean(aux) - param)^2 + se^2,
                      coverage = sum(lo < param & param <= hi)/length(lo)*100)
  row.names(resul) <- ""
  resul
}

print_ptrans <- function(models, from, to, param){
  aux <- unlist(lapply(models, function(.mod) .mod$ptrans[from,to]))
  se <- sqrt(sum((aux - mean(aux))^2)/(length(aux) - 1))
  # lo <- unlist(lapply(models, function(.mod) .mod$qtrans$L[from,to]))
  # hi <- unlist(lapply(models, function(.mod) .mod$qtrans$U[from,to]))
  resul <- data.frame(mean = mean(aux), median = median(aux), 
                      p25 = quantile(aux, probs = 0.025), p975 = quantile(aux, probs = 0.975),
                      bias = mean(aux) - param, bias_tp = (mean(aux) - param)/param*100,
                      accuracy = (mean(aux) - param)^2 + se^2)#,
                      # coverage = sum(lo < param & param <= hi)/length(lo)*100)
  row.names(resul) <- ""
  resul
}

print_HR <- function(models, from = 1, param){
  aux <- unlist(lapply(models, function(.mod) .mod$HZ$FP[from, 'HR']))
  se <- sqrt(sum((aux - mean(aux))^2)/(length(aux) - 1))
  lo <- unlist(lapply(models, function(.mod) .mod$HZ$FP[from, ]['L']))
  hi <- unlist(lapply(models, function(.mod) .mod$HZ$FP[from, ]['U']))
  resul <- data.frame(mean = mean(aux), median = median(aux), 
                      p25 = quantile(aux, probs = 0.025), p975 = quantile(aux, probs = 0.975),
                      bias = mean(aux) - param, bias_tp = (mean(aux) - param)/param*100,
                      accuracy = (mean(aux) - param)^2 + se^2,
                      coverage = sum(lo < param & param <= hi, na.rm = T)/
                        (length(lo) - sum(is.na(lo) | is.na(hi)))*100)
  row.names(resul) <- ""
  resul
}

print_cox <- function(models, tipus, param){
  ind2 <- 'FP'
  if (!is.null(tipus)){
    if (tipus == 1) ind2 <- 'strata(type)scd:FP'
    if (tipus == 2) ind2 <- 'strata(type)ic:FP'
  }
  aux <- unlist(lapply(models, function(.mod) exp(sum(.mod$coef[[ind2]], na.rm = T))))
  se <- sqrt(sum((aux - mean(aux))^2)/(length(aux) - 1))
  lo <- unlist(lapply(models, function(.mod) exp(sum(.mod$confint[ind2, 1], na.rm = T)) ))
  hi <- unlist(lapply(models, function(.mod) exp(sum(.mod$confint[ind2, 2], na.rm = T)) ))
  resul <- data.frame(mean = mean(aux), median = median(aux), 
                      p25 = quantile(aux, probs = 0.025), p975 = quantile(aux, probs = 0.975),
                      bias = mean(aux) - param, bias_tp = (mean(aux) - param)/param*100,
                      accuracy = (mean(aux) - param)^2 + se^2,
                      coverage = sum(lo < param & param <= hi, na.rm = T)/
                        (length(lo) - sum(is.na(lo) | is.na(hi)))*100)
  row.names(resul) <- ""
  resul
}

print_disc <- function(models, ind, param){
  if (!is.null(ind)){
    if (ind == 1) ind2 = 'typescd:FP'
    if (ind == 2) ind2 = 'typeic:FP'
    aux <- unlist(lapply(models, function(.mod) exp(sum(.mod$coef[ind2], na.rm = T))))
    lo <- unlist(lapply(models, function(.mod) exp(sum(.mod$confint[ind2, 1], na.rm = T))))
    hi <- unlist(lapply(models, function(.mod) exp(sum(.mod$confint[ind2, 2], na.rm = T))))
  } else {
    aux <- unlist(lapply(models, function(.mod) exp(sum(.mod$coef['FP'], na.rm = T))))
    lo <- unlist(lapply(models, function(.mod) exp(sum(.mod$confint['FP', 1], na.rm = T)) ))
    hi <- unlist(lapply(models, function(.mod) exp(sum(.mod$confint['FP', 2], na.rm = T)) ))
  }
  se <- sqrt(sum((aux - mean(aux))^2)/(length(aux) - 1))
  resul <- data.frame(mean = mean(aux), median = median(aux), 
                      p25 = quantile(aux, probs = 0.025), p975 = quantile(aux, probs = 0.975),
                      bias = mean(aux) - param, bias_tp = (mean(aux) - param)/param*100,
                      accuracy = (mean(aux) - param)^2 + se^2,
                      coverage = sum(lo < param & param <= hi, na.rm = T)/
                        (length(lo) - sum(is.na(lo) | is.na(hi)))*100)
  row.names(resul) <- ""
  resul
}
