library(dplyr)
set.seed(10); 
n_seed <- sample.int(n = 1000, size = size)

load(sprintf(PATTERN, n_seed[1]))
print(n_seed[1])
msm_50.def <- msm_50

for (.seed in n_seed[2:size]){
  print(.seed)
  load(sprintf(PATTERN, .seed))
  ln <- length(msm_50.def)
  ln2 <- length(msm_50)
  for (.i in 1:ln2){
    msm_50.def[[ln + .i]] <- msm_50[[.i]]
  }
}

msm_50 <- msm_50.def

save(msm_50, file = OFILE)
