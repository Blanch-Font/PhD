#Afegir sensibilitat 85%
library(msm)

#Nsim <- 500 # number of simulation runs
Nsim <- 250 # number of simulation runs
N <- 62000 # number of women in the study
Fup <- 11
left <- F # sense left-truncation

# scale_param <- 688.1962
# shape_param <- 1
exp1_param <- (1/scale_param)^shape_param
exp2_param <- 1/4 #paràmetre exponencial

#msm_50 <- list() # list of MSM models without left-truncation

res_com <- list()
res_obs <- list()

# set.seed(n_seed[.nsim])
for (i in 1:Nsim){
  # Uniform for exponential distribution (Tpre)
  u <- runif(N)
  if (trans12_exp){
    Tpre <- -log(u)/exp1_param
  } else{
    Tpre <- (-log(u)/exp1_param)^(1/shape_param)
  }
  # sojourn time exponential distribution
  u_soj <- runif(N)
  Tsoj <- -log(u_soj)/exp2_param
  Tenter <- 0
  
  # create data set
  # assume fixed screening intervals and every women enters at 50
  my.Data <- data.frame(id = rep(1:N, each = Fup),
                        N_crib = rep(1:Fup, N),
                        Tentry = Tenter, 
                        Tpre = rep(Tpre, each = Fup),
                        Tsoj = rep(Tsoj, each = Fup),
                        obs_type = 1)
  rm(u, u_soj, Tsoj, Tenter)
  
  # time of clinical entrance
  my.Data$Tclin <- with(my.Data, Tpre + Tsoj)
  
  ## time at visits
  # my.Data$time <- with(my.Data, (N_crib - 1)*2) # time at visit
  my.Data$time <- with(my.Data, Tentry*left + (N_crib - 1)*2)
  
  #Eliminem els prevalents
  my.Data <- subset(my.Data, !(Tpre <= Tentry))
  
  ## state variable at visit; 0, 1 or 2: no cancer, pre-clinical cancer, and clinical cancer
  my.Data$state <- with(my.Data, as.numeric(Tpre < time) + as.numeric(Tclin < time))
  
  ## remove screening observations after Tclin+2
  ## hence Tclin < 6 included in data, even if not observed
  ## in order to see how estimates perform if follow-up continued after screen detected cancer
  my.Data <- subset(my.Data, time <= (Tclin+2))
  
  ## by choosing time > Tpre+2, the interval cancers have select=1
  my.Data$observat <- with(my.Data, time <= (Tpre + 2))
  
  ## adjust time in last row if it refers to clinical cancer
  sel_clin <- my.Data$state == 2
  my.Data$time[sel_clin] <- my.Data$Tclin[sel_clin]
  my.Data$obs_type[sel_clin] <- 3
  
  #Cut time at 20 years
  my.Data <- subset(my.Data, time <= 20)
  
  amy.Data <- my.Data
  
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
  
  ##### Initial values transition matrix
  Q <- rbind ( c(-0.004, 0.004, 0), 
               c(0,   -0.55, 0.55),
               c(0,    0,   0) )
  # Q.crude <- crudeinits.msm(state ~ time, subject = id, data = my.Data, qmatrix = Q)
  if (pci){
    if (fixed_pci){
      mod1.com <- msm(formula = state ~ time, subject = id, data = my.Data, qmatrix = Q,
                      obstype = obs_type, pci = c(4, 8, 12, 16), fixedpars = c(4, 6, 8, 10),
                      method = "BFGS", control = list(fnscale = 10000, maxit = 10000))
      # sojourn.msm(mod1.com)
      mod1.obs <- msm(formula = state ~ time, subject = id, data = subset(my.Data, observat),
                      qmatrix = Q, obstype = obs_type, pci = c(4, 8, 12, 16),
                      fixedpars = c(4, 6, 8, 10),
                      method = "BFGS", control = list(fnscale = 10000, maxit = 10000))
      # sojourn.msm(mod1.com)
    } else{
      mod1.com <- msm(formula = state ~ time, subject = id, data = my.Data, qmatrix = Q,
                      obstype = obs_type, pci = c(4, 8, 12, 16),
                      method = "BFGS", control = list(fnscale = 10000, maxit = 10000))
      # sojourn.msm(mod1.com)
      mod1.obs <- msm(formula = state ~ time, subject = id, data = subset(my.Data, observat),
                      qmatrix = Q, obstype = obs_type, pci = c(4, 8, 12, 16),
                      method = "BFGS", control = list(fnscale = 10000, maxit = 10000))
      # sojourn.msm(mod1.com)
    }
  } else{
    mod1.com <- msm(state ~ time, subject = id, data = my.Data, qmatrix = Q, obstype = obs_type,
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