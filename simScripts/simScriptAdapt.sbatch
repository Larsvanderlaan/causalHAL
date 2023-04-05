#!/bin/bash

#SBATCH --array=1
#SBATCH --nodes=2
#SBATCH --output=Rout/par-%J.out
#SBATCH --error=Rout/par-%J.err
#SBATCH --cpus-per-task=1
echo "LOADING R"
module load R
echo "R LOADED"
Rscript -e 'source("R_setup.R"); source("~/causalHAL/simsScriptAdapt.R"); n = as.numeric(Sys.getenv("n")); pos_const = numeric(Sys.getenv("const")); muIsHard = as.logical(Sys.getenv("hard")); do_sims(2000, n, pos_const, muIsHard)'
