#!/bin/bash

#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --output=Rout/par-%J.out
#SBATCH --error=Rout/par-%J.err
echo "LOADING R"
module load R
echo "R LOADED"
Rscript -e 'source("~/causalHAL/scripts/R_setup.R"); source("~/causalHAL/scripts/DR2.R"); n = as.numeric(Sys.getenv("n")); pos_const = as.numeric(Sys.getenv("const")); do_sims(n, pos_const, 5000)'
