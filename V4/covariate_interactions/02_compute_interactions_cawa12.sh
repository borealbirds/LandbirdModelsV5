#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=92G
#SBATCH --time=24:00:00
#SBATCH --job-name=cawa12
#SBATCH --mail-user=mannfred@ualberta.ca

module load StdEnv/2023
module load gcc/9.3.0
module load gdal/3.5.1
module load r/4.4.0

export NODELIST=$(echo $(srun hostname))
Rscript --vanilla 04_compute_interactions_cawa12.R
