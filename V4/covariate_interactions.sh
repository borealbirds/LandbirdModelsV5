#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=92G
#SBATCH --time=01:45:00
#SBATCH --job-name=covariate_interactions_take1
#SBATCH --mail-user=mannfred@ualberta.ca
#SBATCH --mail-type=ALL

module load StdEnv/2023
module load r/4.4.0

export NODELIST=$(echo $(srun hostname))
Rscript --vanilla covariate_interactions_conw_bcr60.R