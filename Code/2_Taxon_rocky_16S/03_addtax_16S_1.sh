#!/bin/sh

#SBATCH --job-name=16S_taxon_1
#SBATCH --output=16S_taxon_1.log
#SBATCH --error=16S_taxon_1.error
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=80
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rwoolive@utk.edu
#SBATCH --time=3-00:00:00

module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2

Rscript ./rocky_dada_16S/03_addtax_16S_1.R \ \
 

