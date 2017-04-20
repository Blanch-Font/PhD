library(dplyr)
set.seed(10); 
n_seed <- sample.int(n = 1000, size = size)
# n_seed <- c(508, 307, 427, 692, 85, 225, 273, 271, 611, 426, 646, 562, 113, 589, 354, 423,
#             52, 260, 821)

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

msm_mod <- msm_mod.def
msm_mod_4 <- msm_mod_4.def
cox_mod <- cox_mod.def
cox_mod_4 <- cox_mod_4.def
disc_mod <- disc_mod.def
disc_mod_4 <- disc_mod_4.def

save(msm_mod, msm_mod_4, cox_mod, cox_mod_4, disc_mod, disc_mod_4, file = OFILE)
