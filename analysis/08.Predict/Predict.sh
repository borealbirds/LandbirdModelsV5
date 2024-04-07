#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=125G
#SBATCH --time=0:01:00
#SBATCH --job-name=NM5_predict
#SBATCH --mail-user=ecknight@ualberta.ca

module load StdEnv/2020
module load gcc/9.3.0
module load gdal/3.5.1
module load r/4.2.1

export NODESLIST=$(echo $(srun hostname))
Rscript --vanille 08.Predict.R
