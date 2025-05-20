
# Original code:
# https://benjjneb.github.io/dada2/assign.html

library(dada2); packageVersion("dada2")
library(ape)





################# Taxonomic assignment: Bacteria ################# 
sampdate <- c("All_agmicrobiome") # fall 2020, spring 2021, summer 2021, fall 2021, spring 2022, summer 2022
region <- "16S"
name <- "Bacteria"



### OPTION 1: assignTaxonomy function


### perform on Rocky because it is computationally intensive:
# connect to UT VPN
# move 4seqtab_16S.csv file to Taxon_rocky_16S folder and copy the folder to the desktop
# in terminal:
# ssh rwoolive@rocky.nimbios.org # if it asks for it, the passphrase is just my macbook password
# on my computer's terminal, copy the folder from my local Desktop to rocky:
# scp -r /Users/rachelwooliver/Desktop/Taxon_rocky_16S rwoolive@rocky.nimbios.org:~/
# start the run: 
# sbatch sbatch Taxon_rocky_16S/03_addtax_16S_1.sh
# check that it's running
# squeue
# if you need to cancel the job:
# scancel [job number]
# if you want to look at the job log:
# less 16S_taxon_1.log
# press q and enter to exit
# if you want to look at the job errors:
# less 16S_taxon_1.error
# output will be in the output folder
# when the job is done, in the terminal for my computer (not rocky), copy from rocky to local desktop
# scp -r rwoolive@rocky.nimbios.org:~/Taxon_rocky_16S/ /Users/rachelwooliver/Desktop/
# copy the output to taxonomy folder
# use control-D to log out of rocky



### OPTION 2
# Use the RDP Classifier
# RDP Naive Bayesian rRNA Classifier Version 2.11, September 2015
# Go to http://rdp.cme.msu.edu/classifier/classifier.jsp
# Citation:  Wang, Q, G. M. Garrity, J. M. Tiedje, and J. R. Cole. 2007. Naïve Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. Appl Environ Microbiol. 73(16):5261-7.

# Under "Choose a gene:" select "16S rRNA training set 18"
# Then click the button that says "Choose file". Navigate to the 
# fasta file for ASVs (4seqtab.fa). Then click "Submit".

# Select confidence threshold of 95%
# click on "show assignment detail for Root only "
# Save the fixed rank output file in the taxonomy folder.

# Open the fix rank file, copy and paste to an excel sheet with one row
# for headers, and save it as a csv.

# Read in the file and trim the dataframe to the asvs that are in the 
# appropriate kingdom.


taxa3 <- readRDS(paste0("Raw-data/sequence/",sampdate,"/",name,"/Taxon_rocky_16S/output/Taxon_rocky_16S_tax_assignation_1.rds"))
#colnames(taxa3) <- gsub("tax.", "", colnames(taxa3))


# remove eukaryotes
table(taxa3$Kingdom) 
taxa3[which(taxa3$Kingdom %in% c("Eukaryota", NA)),] # from corn and unidentifiable
nonbac <- which(taxa3$Kingdom %in% c("Eukaryota", NA))
write.csv(taxa3[nonbac,], paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/2_eukaryotes-removed.csv"))
taxa5 <- taxa3[-nonbac,] # remove the eukaryota and NAs



# remove the same asvs from the seqtab file
seqs <- read.csv(paste0("Raw-data/sequence/",sampdate,"/",name,"/seqtab/4seqtab_16S.csv"))
seqs <- seqs[-nonbac,]
write.csv(seqs, paste0("Raw-data/sequence/",sampdate,"/",name,"/seqtab/5seqtab.csv"))

# remove the same asvs from the fasta file
bacfasta3 <- ape::read.FASTA(paste0("Raw-data/sequence/",sampdate,"/",name,"/seqtab/4seqtab_16S.fa"))
bacfasta5 <- bacfasta3[-nonbac]
ape::write.FASTA(bacfasta5, paste0("Raw-data/sequence/",sampdate,"/",name,"/seqtab/5seqtab.fa"))

# remove the same asvs from the counts file
abund <- read.csv(paste0("Raw-data/sequence/",sampdate,"/",name,"/seqtab/4counts.tsv"), sep="\t", row.names = 1)
abund2 <- abund[-nonbac,]
write.csv(abund2, paste0("Raw-data/sequence/",sampdate,"/",name,"/seqtab/5counts.csv"))



# check that there are no unidentifiable asvs
length(which(is.na(taxa5$Kingdom)==TRUE)) 


# export taxonomy info
write.csv(taxa5, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/3_taxonomy-assignment.csv"))
taxa5 <- read.csv(paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/3_taxonomy-assignment.csv"))


# number of asvs
write.csv(dim(taxa5)[1], paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/3_number-taxa.csv"))




# look at kingdom represented
kingdom <- as.data.frame(table(taxa5$Kingdom))
kingdom$percentage <- round((kingdom$Freq/dim(taxa5)[1])*100,2)
write.csv(kingdom,  paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/4_kingdom.csv"))

# look at phyla represented
phyla <- as.data.frame(table(taxa5$Phylum))
phyla$percentage <- round((phyla$Freq/dim(taxa5)[1])*100,2)
write.csv(phyla,  paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/4a_phyla.csv"))

# look at classes represented
classes <- as.data.frame(table(taxa5$Class))
classes$percentage <- round((classes$Freq/dim(taxa5)[1])*100,2)
write.csv(classes, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/4b_classes.csv"))

# look at orders represented
orders <- as.data.frame(table(taxa5$Order))
orders$percentage <- round((orders$Freq/dim(taxa5)[1])*100,2)
write.csv(orders, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/4c_orders.csv"))

# look at families represented
family <- as.data.frame(table(taxa5$Family))
family$percentage <- round((family$Freq/dim(taxa5)[1])*100,2)
write.csv(family, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/4d_families.csv"))

# look at genera represented
genus <- as.data.frame(table(taxa5$Genus))
genus$percentage <- round((genus$Freq/dim(taxa5)[1])*100,2)
write.csv(genus, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/4e_genera.csv"))


rm(list=ls())















################# Taxonomic assignment: Fungi ################# 
sampdate <- c("All_agmicrobiome") # fall 2020, spring 2021, summer 2021, fall 2021, spring 2022, summer 2022
region <- "ITS"
name <- "Fungi"



### perform on Rocky because it is computationally intensive:
# move 4seqtab_ITS.csv file to Taxon_rocky_ITS folder and copy the folder to the desktop
# in terminal:
# ssh rwoolive@rocky.nimbios.org # the passphrase is just my macbook password
# on my computer's terminal, copy the folder from my local Desktop to rocky:
# scp -r /Users/rachelwooliver/Desktop/Taxon_rocky_ITS rwoolive@rocky.nimbios.org:~/
# start the run: 
# sbatch sbatch Taxon_rocky_ITS/03_addtax_ITS_1.sh
# check that it's running
# squeue
# if you need to cancel the job:
# scancel [job number]
# if you want to look at the job log:
# less ITS_taxon_1.log
# if you want to look at the job errors:
# less ITS_taxon_1.error
# output will be in the output folder
# when the job is done, in the terminal for my computer (not rocky), copy from rocky to local desktop
# scp -r rwoolive@rocky.nimbios.org:~/Taxon_rocky_ITS/ /Users/rachelwooliver/Desktop/
# move the output to All_agmicrobiome/Fungi/taxonomy/
# use control-D to log out of rocky





# OPTION 2:OPTION 2: Use the RDP Classifier
# RDP Naive Bayesian rRNA Classifier Version 2.11, September 2015
# Go to http://rdp.cme.msu.edu/classifier/classifier.jsp
# Citation:  Wang, Q, G. M. Garrity, J. M. Tiedje, and J. R. Cole. 2007. Naïve Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. Appl Environ Microbiol. 73(16):5261-7.

# if 4seqtab.fa is to large, you may need to split into multiple smaller files

# Under "Choose a gene:" select the UNITE database
# Then click the button that says "Choose file". Navigate to the 
# fasta file for ASVs (4seqtab.fa). Then click "Submit".

# Select confidence threshold of 95%
# click on "show assignment detail for Root only "
# Save the fixed rank output file in the taxonomy folder.

# Open the fix rank file, copy and paste to an excel sheet with one row
# for headers, and save it as a csv.

# Read in the file and trim the dataframe to the asvs that are in the 
# appropriate kingdom.



taxa3 <- readRDS(paste0("Raw-data/sequence/",sampdate,"/",name,"/Taxon_rocky_ITS/output/Taxon_rocky_ITS_tax_assignation_1.rds"))
#colnames(taxa3) <- gsub(pattern = "tax.", replacement = "", colnames(taxa3))
table(taxa3$Phylum)

# remove non-fungi
table(taxa3$Kingdom) # no non-fungal ASVs
# write.csv(taxa3, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/2_nonfungi-removed.csv"))
taxa5 <- taxa3#[-nonbac,] # remove the corns



seqs <- read.csv(paste0("Raw-data/sequence/",sampdate,"/",name,"/seqtab/4seqtab_ITS.csv"))
seqs <- seqs#[-nonbac,]
write.csv(seqs, paste0("Raw-data/sequence/",sampdate,"/",name,"/seqtab/5seqtab.csv"))

# remove the same asvs from the fasta file
bacfasta3 <- ape::read.FASTA(paste0("Raw-data/sequence/",sampdate,"/",name,"/seqtab/4seqtab_ITS.fa"))
bacfasta5 <- bacfasta3#[-nonbac]
ape::write.FASTA(bacfasta5, paste0("Raw-data/sequence/",sampdate,"/",name,"/seqtab/5seqtab.fa"))

# remove the same asvs from the counts file
abund <- read.table(paste0("Raw-data/sequence/",sampdate,"/",name,"/seqtab/4counts.tsv") )
abund2 <- abund#[-nonbac,]
write.csv(abund2, paste0("Raw-data/sequence/",sampdate,"/",name,"/seqtab/5counts.csv"))



# check that there are no unidentifiable asvs
length(which(is.na(taxa5$Kingdom)==TRUE)) 


# export taxonomy info
write.csv(taxa5, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/3_taxonomy-assignment.csv"))
taxa5 <- read.csv(paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/3_taxonomy-assignment.csv"))


# number of asvs
write.csv(dim(taxa5)[1], paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/3_number-taxa.csv"))




# look at kingdoms represented
kingdoms <- as.data.frame(table(taxa5$Kingdom))
kingdoms$percentage <- round((kingdoms$Freq/dim(taxa5)[1])*100,2)
write.csv(kingdoms,  paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/4a_kingdoms.csv"))

# look at phyla represented
phyla <- as.data.frame(table(taxa5$Phylum))
phyla$percentage <- round((phyla$Freq/dim(taxa5)[1])*100,2)
write.csv(phyla,  paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/4a_phyla.csv"))

# look at classes represented
classes <- as.data.frame(table(taxa5$Class))
classes$percentage <- round((classes$Freq/dim(taxa5)[1])*100,2)
write.csv(classes, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/4b_classes.csv"))

# look at orders represented
orders <- as.data.frame(table(taxa5$Order))
orders$percentage <- round((orders$Freq/dim(taxa5)[1])*100,2)
write.csv(orders, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/4c_orders.csv"))

# look at families represented
family <- as.data.frame(table(taxa5$Family))
family$percentage <- round((family$Freq/dim(taxa5)[1])*100,2)
write.csv(family, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/4d_families.csv"))

# look at genera represented
genus <- as.data.frame(table(taxa5$Genus))
genus$percentage <- round((genus$Freq/dim(taxa5)[1])*100,2)
write.csv(genus, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/4e_genera.csv"))

# look at species represented
taxa5$Genus <- gsub("g__", "", taxa5$Genus)
taxa5$Species <- gsub("s__", "", taxa5$Species)
species <- as.data.frame(table(paste(taxa5$Genus, taxa5$Species)))
species$percentage <- round((species$Freq/dim(taxa5)[1])*100,2)
write.csv(species, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/4e_species.csv"))



########### Assign function using Functional trait database

ftrait <- read.csv("Raw-data/Fungal-trait-database/data 2.csv")

taxa <- read.csv(paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/3_taxonomy-assignment.csv"))
taxa$Genus <- gsub(pattern="g__", replacement="", taxa$Genus)
taxa$Primary_lifestyle <- rep(NA, dim(taxa)[1])
taxa$Secondary_lifestyle <- rep(NA, dim(taxa)[1])

for(i in 1:dim(taxa)[1]){
  if(taxa$Genus[i] %in% ftrait$GENUS){
    taxa$Primary_lifestyle[i] <- ftrait$primary_lifestyle[which(ftrait$GENUS==taxa$Genus[i])]
    taxa$Secondary_lifestyle[i] <- ftrait$Secondary_lifestyle[which(ftrait$GENUS==taxa$Genus[i])]
  }
}


# functional trait distribution
ftraittab <- as.data.frame(table(taxa$Primary_lifestyle))
ftraittab$percentage <- round((ftraittab$Freq/dim(taxa5)[1])*100,2)
write.csv(ftraittab, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/5a_primary-lifestyle-asvs.csv"))

# create data subsets for certain functional groups
write.csv(taxa[which(taxa$Primary_lifestyle=="arbuscular_mycorrhizal"),], paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/5b_amf-asvs.csv"))
write.csv(dim(taxa[which(taxa$Primary_lifestyle=="arbuscular_mycorrhizal"),])[1], paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/5b_amf-asvs-num.csv"))
write.csv(taxa[which(taxa$Primary_lifestyle=="soil_saprotroph"),], paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/5b_soil_saprotroph-asvs.csv"))
write.csv(dim(taxa[which(taxa$Primary_lifestyle=="soil_saprotroph"),])[1], paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/5b_soil_saprotroph-asvs-num.csv"))
write.csv(taxa[which(taxa$Primary_lifestyle=="litter_saprotroph"),], paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/5b_litter_saprotroph-asvs.csv"))
write.csv(dim(taxa[which(taxa$Primary_lifestyle=="litter_saprotroph"),])[1], paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/5b_litter_saprotroph-asvs-num.csv"))
write.csv(taxa[which(taxa$Primary_lifestyle=="plant_pathogen"),], paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/5b_plant-pathogen-asvs.csv"))
write.csv(dim(taxa[which(taxa$Primary_lifestyle=="plant_pathogen"),])[1], paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/5b_plant_pathogen-asvs-num.csv"))


write.csv(taxa, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/5_taxonomy-assignment-with-ft.csv"))




#### Subset to glomeromycota to explore AMF


taxa5 <- read.csv(paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/5_taxonomy-assignment-with-ft.csv"))
amftax <- taxa5[which(taxa5$Phylum=="p__Glomeromycota"),]

write.csv(amftax, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/6_amf-Glomeromycota.csv"))
write.csv(dim(amftax)[1], paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/6_amf-Glomeromycota-num.csv"))




# look at phyla represented
phyla <- as.data.frame(table(amftax$Phylum))
phyla$percentage <- round((phyla$Freq/dim(amftax)[1])*100,2)
write.csv(phyla,  paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/6a_amfphyla.csv"))

# look at classes represented
classes <- as.data.frame(table(amftax$Class))
classes$percentage <- round((classes$Freq/dim(amftax)[1])*100,2)
write.csv(classes, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/6b_amfclasses.csv"))

# look at orders represented
orders <- as.data.frame(table(amftax$Order))
orders$percentage <- round((orders$Freq/dim(amftax)[1])*100,2)
write.csv(orders, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/6c_amforders.csv"))

# look at families represented
family <- as.data.frame(table(amftax$Family))
family$percentage <- round((family$Freq/dim(amftax)[1])*100,2)
write.csv(family, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/6d_amffamilies.csv"))

# look at genera represented
genus <- as.data.frame(table(amftax$Genus))
genus$percentage <- round((genus$Freq/dim(amftax)[1])*100,2)
write.csv(genus, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/6e_amfgenera.csv"))

# look at species represented
amftax$Species <- gsub(pattern="s__", replacement="", amftax$Species)
species <- as.data.frame(table(paste(amftax$Genus, amftax$Species)))
species$percentage <- round((species$Freq/dim(amftax)[1])*100,2)
write.csv(species, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/6e_amfspecies.csv"))





#### Explore plant pathogenic fungi


taxa5 <- read.csv(paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/5_taxonomy-assignment-with-ft.csv"))
pathtax <- taxa5[which(taxa5$Primary_lifestyle=="plant_pathogen"),]


# look at phyla represented
phyla <- as.data.frame(table(pathtax$Phylum))
phyla$percentage <- round((phyla$Freq/dim(pathtax)[1])*100,2)
write.csv(phyla,  paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/5c_plant-pathogen-phyla.csv"))

# look at classes represented
classes <- as.data.frame(table(pathtax$Class))
classes$percentage <- round((classes$Freq/dim(pathtax)[1])*100,2)
write.csv(classes, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/5c_plant-pathogen-classes.csv"))

# look at orders represented
orders <- as.data.frame(table(pathtax$Order))
orders$percentage <- round((orders$Freq/dim(pathtax)[1])*100,2)
write.csv(orders, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/5c_plant-pathogen-orders.csv"))

# look at families represented
family <- as.data.frame(table(pathtax$Family))
family$percentage <- round((family$Freq/dim(pathtax)[1])*100,2)
write.csv(family, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/5c_plant-pathogen-families.csv"))

# look at genera represented
genus <- as.data.frame(table(pathtax$Genus))
genus$percentage <- round((genus$Freq/dim(pathtax)[1])*100,2)
write.csv(genus, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/5c_plant-pathogen-genera.csv"))

# look at species represented
pathtax$Species <- gsub(pattern="s__", replacement="", pathtax$Species)
species <- as.data.frame(table(paste(pathtax$Genus, pathtax$Species)))
species$percentage <- round((species$Freq/dim(pathtax)[1])*100,2)
write.csv(species, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/5c_plant-pathogen-species.csv"))



