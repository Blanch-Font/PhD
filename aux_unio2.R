library(dplyr)
set.seed(10); 
n_seed <- sample.int(n = 1000, size = size)[ini:size]

load(sprintf(PATTERN, n_seed[1]))
print(n_seed[1])
msm_50.def <- msm_50
msm_left.def <- msm_left
cox_50.def <- cox_50
cox_left.def <- cox_left
msm_50_4.def <- msm_50_4
msm_left_4.def <- msm_left_4
cox_50_4.def <- cox_50_4
cox_left_4.def <- cox_left_4
disc_50.def <- disc_50
disc_left.def <- disc_left
disc_50_4.def <- disc_50_4
disc_left_4.def <- disc_left_4

for (.seed in n_seed[2:(size - (ini - 1))]){
  print(.seed)
  load(sprintf(PATTERN, .seed))
  ln <- length(msm_50.def)
  ln2 <- length(msm_50)
  for (.i in 1:ln2){
    msm_50.def[[ln + .i]] <- msm_50[[.i]]
    msm_left.def[[ln + .i]] <- msm_left[[.i]]
    cox_50.def[[ln + .i]] <- cox_50[[.i]]
    cox_left.def[[ln + .i]] <- cox_left[[.i]]
    msm_50_4.def[[ln + .i]] <- msm_50_4[[.i]]
    msm_left_4.def[[ln + .i]] <- msm_left_4[[.i]]
    cox_50_4.def[[ln + .i]] <- cox_50_4[[.i]]
    cox_left_4.def[[ln + .i]] <- cox_left_4[[.i]]
    disc_50.def[[ln + .i]] <- disc_50[[.i]]
    disc_left.def[[ln + .i]] <- disc_left[[.i]]
    disc_50_4.def[[ln + .i]] <- disc_50_4[[.i]]
    disc_left_4.def[[ln + .i]] <- disc_left_4[[.i]]
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
