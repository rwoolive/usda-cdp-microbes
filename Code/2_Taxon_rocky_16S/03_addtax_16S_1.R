## ------------------------------------------------------------------------
#This script adds taxonomy to each of the inferred ASVs

library(dada2)
library(tidyverse)


cat(paste0('\n',"You are using DADA2 version ", packageVersion('dada2'),'\n'))

cat('################################\n\n')

args <- commandArgs(trailingOnly = TRUE)

seqtab.nochim <- "Taxon_rocky_16S/4seqtab_16S.csv" #change to your file structure
output_16S <- "Taxon_rocky_16S/output/"
name <- "Taxon_rocky_16S"
tax_db <- "Taxon_rocky_16S/rdp_train_set_18.fa.gz" #change to your tax database (Silva, UNITE, MaarjAM*)
method <- "dada"
threshold <- as.integer(50)

dir.create(file.path(output_16S), showWarnings = FALSE)
# dir.create(file.path(output_16S, "03_taxonomy"), showWarnings = FALSE)
# dir.create(file.path(output_16S, "03_taxonomy", name), showWarnings = FALSE)

#output_16S <- paste0(output_16S,"/03_taxonomy/",name,"/")

# Assign taxonomy (general)
set.seed(42) #random  generator necessary for reproducibility

seqtab.nochim <- read.csv(seqtab.nochim)

#seqtab.nochim <- seqtab.nochim[,1:200000]
# seqtab.nochim <- seqtab.nochim[,1:5] #test


# single step
#taxid <- assignTaxonomy(colnames(seqtab.nochim)[1], tax_db, multithread = TRUE, tryRC = TRUE)

# multi step


num_seq = dim(seqtab.nochim)
num_step = 10
seq_per_step = as.integer(num_seq[1] / num_step)

cat('\n')
cat(paste0('# Setting seqs per step to  "', seq_per_step, '"\n'))
cat('\n# And beginning taxon assignment. \n\n')


last_step_num_seq = num_seq[1] %% num_step

taxa_0 = data.frame(Kingdom=NA,Phylum=NA,Class=NA,Order=NA,Family=NA,Genus=NA)#,Species=NA)
for(i in 1:num_step){
  sbegin = 1 + (i-1)*seq_per_step
  send = i*seq_per_step
  print(paste("Processing sequences", sbegin, "through", send, "of", num_seq[1]))
  taxa_s <- assignTaxonomy(seqtab.nochim$X[sbegin:send], tax_db, multithread=TRUE, tryRC = TRUE)
  taxa_0 <-rbind(taxa_0, as.data.frame(taxa_s))
  print(paste("Finished Step", i, "of", num_step))
}
sbegin = 1 + num_step*seq_per_step
send = num_seq[1]
print(paste("Processing sequences", sbegin, "through", send, "of", num_seq[1]))
taxa_s <- assignTaxonomy(seqtab.nochim$X[sbegin:send], tax_db, multithread=TRUE, tryRC = TRUE)
taxa_0 <-rbind(taxa_0, as.data.frame(taxa_s))
print("Finished Remainder")

taxa_0 <- taxa_0[2:dim(taxa_0)[1],]






# Write to disk
saveRDS(taxa_0, paste0(output_16S, name, "_tax_assignation_1.rds"))

cat('\n')
cat(paste0('# The obtained taxonomy file can be found in "', paste0(output_16S, name, "_tax_assignation_1.rds"), '"\n'))
cat('\n# All done!\n\n')
