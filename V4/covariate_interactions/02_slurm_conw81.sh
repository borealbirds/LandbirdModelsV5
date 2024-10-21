#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=16
#SBATCH --ntasks=32
#SBATCH --mem=186G
#SBATCH --time=06:00:00
#SBATCH --job-name=conw81
#SBATCH --mail-user=mannfred@ualberta.ca

module load StdEnv/2023
module load gcc/9.3.0
module load gdal/3.5.1
module load r/4.4.0

export NODELIST=$(echo $(srun hostname))
Rscript --vanilla 04_compute_interactions_conw81.R
