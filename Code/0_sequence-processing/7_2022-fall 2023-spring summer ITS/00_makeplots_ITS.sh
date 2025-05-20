#!/bin/bash


#SBATCH --job-name=rocky_dada_ITS
#SBATCH --output=rocky_dada_ITS.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40

module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2

Rscript ./rocky_dada_ITS/00_makeplots_ITS.R 
