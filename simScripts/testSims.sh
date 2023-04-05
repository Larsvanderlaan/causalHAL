#!/bin/usr/env bash
filename="simsCATE_xgboost"
nsims=2500
export R_LIBS=~/Rlibs2
export R_LIBS_USER=~/Rlibs2
for n in 1000
do
  for const in 3
  do
    sbatch  --export=n=$n,const=$const ~/sieveSims/testSims.sbatch
  done
done
