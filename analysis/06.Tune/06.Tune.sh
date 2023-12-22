#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=0
#SBATCH --time=1:00:00
#SBATCH --job-name=NM5_tuning
#SBATCH --mail-user=ecknight@ualberta.ca

module load StdEnv/2020
module load gcc/9.3.0
module load r/4.2.1

export NODESLIST=$(echo $(srun hostname))
Rscript --vanilla 02.Tune.R
