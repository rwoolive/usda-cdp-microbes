#!/bin/bash


#SBATCH --job-name=rocky_dada_16S
#SBATCH --output=rocky_dada_16S/rocky_dada_16S.log
#SBATCH --error=rocky_dada_16S/rocky_dada_16S.error
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=80
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rwoolive@utk.edu
#SBATCH --time=3-00:00:00




module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2

Rscript ./rocky_dada_16S/02_chiremmer_16S.R \
