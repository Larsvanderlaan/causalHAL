#!/bin/usr/env bash
nsims=2500
export R_LIBS=~/Rlibs2
export R_LIBS_USER=~/Rlibs2
for n in 500 1000 2000 3000 4000 5000
do
  for const in 0.5 1 2 3
  do
    for hard in "TRUE" "FALSE"
    do
    sbatch  --export=n=$n,const=$const,hard=$hard ~/causalHAL/simulationScripts/simScriptAdapt.sbatch
    done
  done
done
