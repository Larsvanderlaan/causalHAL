#!/bin/usr/env bash
filename="simsCATE_xgboost"
nsims=2500
export R_LIBS=~/Rlibs2
export R_LIBS_USER=~/Rlibs2
for n in 500 1000 2000 3000 4000 5000
do
  for const in 3 5 8
  do
    sbatch  --export=n=$n,const=$const ~/sieveSims/simsComplex.sbatch
  done
done
