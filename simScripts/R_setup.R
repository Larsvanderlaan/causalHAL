.libPaths( c( "~/Rlibs2", .libPaths()) )
print(.libPaths())
setwd("~/sieveSims")
print(getwd())
nsims = 5000

library(autocausalML)
#library(future)
#plan(multisession, workers = 16)

