#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=16  
#SBATCH --mem=16G
#SBATCH --time=01:00:00  
#SBATCH --job-name=covariate_interactions_take1
#SBATCH --mail-user=mannfred@ualberta.ca
#SBATCH --mail-type=ALL

module load StdEnd/2020
module load r/4.4.0

export NODELIST=$(echo $(srun hostname))
Rscript --vanilla covariate_interactions.R