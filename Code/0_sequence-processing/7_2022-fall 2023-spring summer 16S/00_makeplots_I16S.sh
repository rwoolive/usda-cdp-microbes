#!/bin/bash


#SBATCH --job-name=rocky_dada_16S
#SBATCH --output=rocky_dada_16S.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40

module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2

Rscript ./rocky_dada_16S/00_makeplots_16S.R
