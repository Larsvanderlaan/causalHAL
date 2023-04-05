#!/bin/bash

#SBATCH --array=1
#SBATCH --nodes=11
#SBATCH --output=Rout/par-%J.out
#SBATCH --error=Rout/par-%J.err
#SBATCH --cpus-per-task=1
echo "LOADING R"
module load R
echo "R LOADED"
Rscript -e 'source("R_setup.R"); n = Sys.getenv("n"); const = Sys.getenv("const"); source("~/sieveSims/simsComplex.R")'
