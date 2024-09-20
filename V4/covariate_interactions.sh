#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=186G
#SBATCH --time=08:00:00
#SBATCH --tmp=100G
#SBATCH --job-name=covariate_interactions_take1
#SBATCH --mail-user=mannfred@ualberta.ca
#SBATCH --mail-type=ALL

module load StdEnv/2023
module load r/4.4.0

rm ~/.RData

export NODELIST=$(echo $(srun hostname))
Rscript --vanilla covariate_interactions.R
