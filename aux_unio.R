library(dplyr)
set.seed(10); 
n_seed <- sample.int(n = 1000, size = size)

load(sprintf(PATTERN, n_seed[1]))
print(n_seed[1])
msm_mod.def <- msm_mod
msm_mod_4.def <- msm_mod_4
cox_mod.def <- cox_mod
cox_mod_4.def <- cox_mod_4
disc_mod.def <- disc_mod
disc_mod_4.def <- disc_mod_4

for (.seed in n_seed[2:size]){
  print(.seed)
  load(sprintf(PATTERN, .seed))
  ln <- length(msm_mod.def)
  ln2 <- length(msm_mod)
  for (.i in 1:ln2){
    msm_mod.def[[ln + .i]] <- msm_mod[[.i]]
    msm_mod_4.def[[ln + .i]] <- msm_mod_4[[.i]]
    cox_mod.def[[ln + .i]] <- cox_mod[[.i]]
    cox_mod_4.def[[ln + .i]] <- cox_mod_4[[.i]]
    disc_mod.def[[ln + .i]] <- disc_mod[[.i]]
    disc_mod_4.def[[ln + .i]] <- disc_mod_4[[.i]]
  }
}

msm_50 <- msm_50.def
msm_left <- msm_left.def
cox_50 <- cox_50.def
cox_left <- cox_left.def
msm_50_4 <- msm_50_4.def
msm_left_4 <- msm_left_4.def
cox_50_4 <- cox_50_4.def
cox_left_4 <- cox_left_4.def
disc_50 <- disc_50.def
disc_left <- disc_left.def
disc_50_4 <- disc_50_4.def
disc_left_4 <- disc_left_4.def

save(msm_50, msm_left, cox_50, cox_left,
     msm_50_4, msm_left_4, cox_50_4, cox_left_4,
     disc_50, disc_left, disc_50_4, disc_left_4, file = OFILE)
