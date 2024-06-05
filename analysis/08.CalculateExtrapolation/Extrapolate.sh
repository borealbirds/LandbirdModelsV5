#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=125G
#SBATCH --time=0:10:00
#SBATCH --job-name=NM5_exptrapolate
#SBATCH --mail-user=ecknight@ualberta.ca

module load StdEnv/2020 gcc/9.3.0 gdal/3.5.1 r/4.2.1 udunits/2.2.28

export NODESLIST=$(echo $(srun hostname))
Rscript --vanilla 08.CalculateExtrapolation.R
