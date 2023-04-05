.libPaths( c( "~/Rlibs2", .libPaths()) )
print(.libPaths())
setwd("~/sieveSims")
print(getwd())
nsims = 1000

library(causalHAL)
#library(future)
#plan(multisession, workers = 16)

