#!/bin/sh

#SBATCH --job-name=ITS_taxon_1
#SBATCH --output=ITS_taxon_1.log
#SBATCH --error=ITS_taxon_1.error
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=80
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rwoolive@utk.edu
#SBATCH --time=3-00:00:00

module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2

Rscript ./Taxon_rocky_ITS/03_addtax_ITS_1.R \ \
 

