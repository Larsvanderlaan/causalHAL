#!/bin/usr/env bash
nsims=2500
export R_LIBS=~/Rlibs2
export R_LIBS_USER=~/Rlibs2
for n in 500 1000 2000 3000 4000 5000
do
  for const in 1 2 3
    for hard in true false
  do
    sbatch  --export=n=$n,const=$const, hard=$hard, ~/causalHAL/simScriptAdapt.sbatch
  done
done
