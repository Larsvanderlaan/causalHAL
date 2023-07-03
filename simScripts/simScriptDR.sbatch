#!/bin/bash

#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=11
#SBATCH --output=Rout/par-%J.out
#SBATCH --error=Rout/par-%J.err
echo "LOADING R"
module load R
echo "R LOADED"
Rscript -e 'source("~/causalHAL/simScripts/R_setup.R"); source("~/causalHAL/simScripts/DR.R"); n = as.numeric(Sys.getenv("n")); pos_const = as.numeric(Sys.getenv("const")); do_sims(n, pos_const, 10000)'
