#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=0
#SBATCH --time=24:00:00
#SBATCH --job-name=NM5_tune_2
#SBATCH --mail-user=ecknight@ualberta.ca
#SBATCH --mail-type=ALL

module load StdEnv/2020 gcc/9.3.0 gdal/3.5.1 r/4.2.1

Rscript 05.Tune.R
