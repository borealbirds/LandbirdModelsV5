#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=2
#SBATCH --ntasks=32
#SBATCH --mem=125G
#SBATCH --time=48:00:00
#SBATCH --job-name=NM5_simplifying
#SBATCH --mail-user=ecknight@ualberta.ca

module load StdEnv/2020
module load gcc/9.3.0
module load gdal/3.5.1
module load r/4.2.1

export NODESLIST=$(echo $(srun hostname))
Rscript --vanilla 07.Simplify.R
