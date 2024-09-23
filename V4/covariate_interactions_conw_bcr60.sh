#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=92G
#SBATCH --time=01:00:00
#SBATCH --job-name=covariate_interactions_take1
#SBATCH --mail-user=mannfred@ualberta.ca

module load StdEnv/2023
module load r/4.4.0

export NODELIST=$(echo $(srun hostname))
Rscript --vanilla covariate_interactions_conw_bcr60.R