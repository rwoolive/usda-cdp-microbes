



##############################################
# DADA2 is a pipeline that processes raw 16S and ITS sequence data. 
# It works entirely in R.
# We have filtered & trimmed, dereplicated, implemented sample inference to 
# determine unique sequences, merged paired reads, and constructed a sequence
# table for each sampling period in each target group.

# In this code, we will merge 16S the sequence tables across sampling periods, 
# then do the same for ITS sequence tables.

# the main tutorial is here: https://benjjneb.github.io/dada2/index.html
# and the tutorial curated for analysis UT genomics core data is here: file:///Users/rachelwooliver/Library/Containers/com.apple.mail/Data/Library/Mail%20Downloads/B273FB41-2BC6-4285-8301-2C3265275660/DADA2Script_2020_12_14.html


##############################################
##############################################
#install DADA2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")
#Bioconductor version 3.12 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

# load libraries
library(dada2)
packageVersion("dada2") # [1] ‘1.30.0’





###################
# Bacteria
###################

### need to replace these samples with redos:
### need to replace these samples with redos:
to_replace <- c("110-16S_2021-07", 
                "118-16S_2021-07", 
                "202-16S_2021-07",
                "215-16S_2021-07",
                "216-16S_2021-07",
                "318-16S_2021-07",
                "219-16S_2021-09",               
                "220-16S_2021-09")
replacement <- c("110r-16S_2021-07", 
                 "118r-16S_2021-07", 
                 "202r-16S_2021-07",
                 "215r-16S_2021-07",
                 "216r-16S_2021-07",
                 "318r-16S_2021-07",
                 "219_16S_2021-fall",
                 "220_16S_2021-fall")

# sample dates
sampdate <- c("2020-09", # fall 2020, 
              "2021-04", # spring 2021, 
              "2021-07", # summer 2021
              "2021-09", # fall 2021
              "2022-spring", # spring 2022
              "2022-summer", # summer 2022
              "2022-fall 2023-spring summer", # fall 2022, spring 2023, some of summer 2023
              "2023 summer-2024 fall") # rest of summer 2023, fall 2023, spring 2024, summer 2024, fall 2024
region <- "16S"
name <- "Bacteria"



# read in data from each sampling period (post-chimera-removal)
tab1 <- t(read.csv(paste0("Raw-data/sequence/",sampdate[1],"_agmicrobiome/",name,"/seqtab/2seqtab.csv"), row.names = 1)) 
tab2 <- t(read.csv(paste0("Raw-data/sequence/",sampdate[2],"_agmicrobiome/",name,"/seqtab/2seqtab.csv"), row.names = 1))
tab3 <- t(read.csv(paste0("Raw-data/sequence/",sampdate[3],"_agmicrobiome/",name,"/seqtab/2seqtab.csv"), row.names = 1))
tab4 <- t(read.csv(paste0("Raw-data/sequence/",sampdate[4],"_agmicrobiome/",name,"/seqtab/2seqtab.csv"), row.names = 1))
tab5 <- t(read.csv(paste0("Raw-data/sequence/",sampdate[5],"_agmicrobiome/",name,"/seqtab/2seqtab.csv"), row.names = 1))
tab6 <- t(read.csv(paste0("Raw-data/sequence/",sampdate[6],"_agmicrobiome/",name,"/seqtab/2seqtab.csv"), row.names = 1))
tab7 <- readRDS(file=paste0("Raw-data/sequence/2022-fall 2023-spring summer_agmicrobiome/Bacteria/rocky_dada_16S/output/02_nochimera_mergeruns/rocky_dada_16S/rocky_dada_16S_seqtab_final_renamed.rds"))
tab8 <- readRDS(file=paste0("Raw-data/sequence/2023 summer-2024 fall_agmicrobiome/16S/rocky_dada_16S/output/02_nochimera_mergeruns/rocky_dada_16S/rocky_dada_16S_seqtab_final_renamed.rds"))
# The sequence table is a matrix with rows corresponding to (and named by) the 
# samples, and columns corresponding to (and named by) the sequence variants. 



# clean up sample names
rownames(tab1) <- gsub(".1.16S", "-16S", rownames(tab1))
rownames(tab1) <- gsub("X", "", rownames(tab1))
rownames(tab1) <- gsub("2020.", "2020-", rownames(tab1))
rownames(tab2) <- gsub(".0.10.16S", "-16S", rownames(tab2))
rownames(tab2) <- gsub("X", "", rownames(tab2))
rownames(tab2) <- gsub("2021.", "2021-", rownames(tab2))
rownames(tab3) <- gsub("_16S", "-16S", rownames(tab3))
rownames(tab3) <- gsub("X", "", rownames(tab3))
rownames(tab3) <- gsub("2021.", "2021-", rownames(tab3))
rownames(tab4) <- gsub("_16S", "-16S", rownames(tab4))
rownames(tab4) <- gsub("X", "", rownames(tab4))
rownames(tab4) <- gsub("2021.", "2021-", rownames(tab4))
rownames(tab5) <- gsub("_16S", "-16S", rownames(tab5))
rownames(tab5) <- gsub("X", "", rownames(tab5))
rownames(tab5) <- gsub("2022.", "2022-", rownames(tab5))
rownames(tab6) <- gsub("_16S", "-16S", rownames(tab6))
rownames(tab6) <- gsub("X", "", rownames(tab6))
rownames(tab6) <- gsub("2022.", "2022-", rownames(tab6))

rownames(tab7) <- gsub("-1", "_16S", rownames(tab7))
rownames(tab7) <- gsub("Fall_2021", "2021-fall", rownames(tab7))
rownames(tab7) <- gsub("Fall_2022", "2022-fall", rownames(tab7))
rownames(tab7) <- gsub("Spring_2023", "2023-spring", rownames(tab7))
rownames(tab7) <- gsub("Summer_2023", "2023-summer", rownames(tab7))

rownames(tab8) <- gsub("-1", "_16S", rownames(tab8))
rownames(tab8) <- gsub("Fall_2021", "2021-fall", rownames(tab8))
rownames(tab8) <- gsub("Fall_2022", "2022-fall", rownames(tab8))
rownames(tab8) <- gsub("Spring_2023", "2023-spring", rownames(tab8))
rownames(tab8) <- gsub("Summer_2023", "2023-summer", rownames(tab8))




# inspect sequence lengths:
table(nchar(colnames(tab1)))
table(nchar(colnames(tab2)))
table(nchar(colnames(tab3)))
table(nchar(colnames(tab4)))
table(nchar(colnames(tab5)))
table(nchar(colnames(tab6)))
table(nchar(colnames(tab7)))
table(nchar(colnames(tab8)))
# Make sure the lengths of the merged sequences fall within the
# expected range for your targeted amplicon.

# NOTE: Sequences that are much longer or shorter than expected may be the 
# result of non-specific priming. You can remove non-target-length sequences 
# from your sequence table: 
# eg. seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256] 
# This is analogous to “cutting a band” in-silico to get amplicons of the 
# targeted length.

# make sure the plots are in order
rownames(tab1)
rownames(tab2)
rownames(tab3)
rownames(tab4)
rownames(tab5)
rownames(tab6)
rownames(tab7)

tab1 <-  tab1[order(rownames(tab1)),]
tab2 <-  tab2[order(rownames(tab2)),]
tab3 <-  tab3[order(rownames(tab3)),]
tab4 <-  tab4[order(rownames(tab4)),]
tab5 <-  tab5[order(rownames(tab5)),]
tab6 <-  tab6[order(rownames(tab6)),]
#tab7 <-  tab7[order(rownames(tab7)),]
#tab8 <-  tab8[order(rownames(tab8)),]

# Combine datasets from different sampling periods:
seqtab_combo <- mergeSequenceTables(tab1, tab2, tab3, tab4, tab5, tab6, as.matrix(tab7), as.matrix(tab8))

# remove bad samples
seqtab_combo_new <- seqtab_combo
# check that we are removing the right samples
dat <- data.frame(a = rownames(seqtab_combo_new[which(rownames(seqtab_combo_new) %in% to_replace),]), 
                  b = to_replace, 
                  c= rownames(seqtab_combo_new[which(rownames(seqtab_combo_new) %in% replacement),]), 
                  d = replacement); dat
# replace
seqtab_combo_new[which(rownames(seqtab_combo_new) %in% to_replace),] <- seqtab_combo_new[which(rownames(seqtab_combo_new) %in% replacement),]
seqtab_combo_new[which(rownames(seqtab_combo_new) %in% c(to_replace[1], replacement[1])),1:2] # check
# remove old columns for replacements
seqtab_combo_new <- seqtab_combo_new[-which(rownames(seqtab_combo_new) %in% replacement),]

seqtab_combo <- seqtab_combo_new

rownames(seqtab_combo)

##############################################
# Remove sequences that are found within the PCR blanks,
# and then remove sequences with less than 10 reads

seqtab.nochim_df <- as.data.frame(t(seqtab_combo))
# determine if any asvs were detected in pcr blanks
unique(seqtab.nochim_df$`PCRblank.16S_2020-09`); which(seqtab.nochim_df$`PCRblank.16S_2020-09`>0) 
unique(seqtab.nochim_df$`Kivlinblank-16S_2021-07`); which(seqtab.nochim_df$`Kivlinblank-16S_2021-07`>0) 
unique(seqtab.nochim_df$`Kivlin_BacteriaBlank-16S_2021-09`); which(seqtab.nochim_df$`Kivlin_BacteriaBlank-16S_2021-09`>0) 
unique(seqtab.nochim_df$`Blank-16S_2022-spring`); which(seqtab.nochim_df$`Blank-16S_2022-spring`>0) 
unique(seqtab.nochim_df$`Blank-16S_2022-summer`); which(seqtab.nochim_df$`Blank-16S_2022-summer`>0) 
unique(seqtab.nochim_df$blank1_Fall_2024); which(seqtab.nochim_df$blank1_Fall_2024>0) 


# there were some asvs detected in the blanks, but they had <10 reads so they  
# will be removed anyway

# remove column with pcr blank
blankcol <- which(colnames(seqtab.nochim_df) %in% c("PCRblank.16S_2020-09", "Kivlinblank-16S_2021-07", "Kivlin_BacteriaBlank-16S_2021-09", "Blank-16S_2022-spring", "Blank-16S_2022-summer", "Undetermined_S0_L001", "blank1_Fall_2024"))
seqtab.nochim_df2 <- seqtab.nochim_df[,-blankcol]


#### remove asvs with less than 10 reads across samples
seqtab.nochim_df2$sum <- rowSums(seqtab.nochim_df2)
length(seqtab.nochim_df2$sum) # number of asvs


# number of asvs with less than 10 reads
write.csv(length(which(seqtab.nochim_df2$sum<10)), file=paste0("Raw-data/sequence/All_agmicrobiome/",name,"/seqtab/3lessthan10reads.csv"))
# number of asvs with > 10 reads
write.csv(length(which(seqtab.nochim_df2$sum>=10)), file=paste0("Raw-data/sequence/All_agmicrobiome/",name,"/seqtab/3greaterthan10reads.csv"))
# proportion of asvs retained
write.csv(length(which(seqtab.nochim_df2$sum>=10))/length(seqtab.nochim_df2$sum), file=paste0("Raw-data/sequence/All_agmicrobiome/",name,"/seqtab/3proportionretained.csv"))
 
# check sample names
colnames(seqtab.nochim_df2)[1:80+80*0]
colnames(seqtab.nochim_df2)[1:80+80*1]
colnames(seqtab.nochim_df2)[1:80+80*2]
colnames(seqtab.nochim_df2)[1:80+80*3]
colnames(seqtab.nochim_df2)[1:80+80*4]
colnames(seqtab.nochim_df2)[1:80+80*5]
colnames(seqtab.nochim_df2)[1:80+80*6]
colnames(seqtab.nochim_df2)[1:80+80*7]
colnames(seqtab.nochim_df2)[1:80+80*8]
colnames(seqtab.nochim_df2)[1:80+80*9]


# remove
seqtab.nochim_df3 <- seqtab.nochim_df2[-which(seqtab.nochim_df2$sum<10),]
hist(seqtab.nochim_df3$sum); range(seqtab.nochim_df3$sum) # check

write.csv(seqtab.nochim_df3, file=paste0("Raw-data/sequence/All_agmicrobiome/",name,"/seqtab/4seqtab_16S.csv"))
seqtab.nochim_df3 <- read.csv(file=paste0("Raw-data/sequence/All_agmicrobiome/Bacteria/seqtab/4seqtab_16S.csv"))


# check sample names
seqtab.nochim_df_check <- read.csv(file=paste0("Raw-data/sequence/All_agmicrobiome/Bacteria/seqtab/4seqtab_16S.csv"))
colnames(seqtab.nochim_df_check)[1+1:80+80*0]
colnames(seqtab.nochim_df_check)[1+1:80+80*1]
colnames(seqtab.nochim_df_check)[1+1:80+80*2]
colnames(seqtab.nochim_df_check)[1+1:80+80*3]
colnames(seqtab.nochim_df_check)[1+1:80+80*4]
colnames(seqtab.nochim_df_check)[1+1:80+80*5]
colnames(seqtab.nochim_df_check)[1+1:80+80*6]
colnames(seqtab.nochim_df_check)[1+1:80+80*7]
colnames(seqtab.nochim_df_check)[1+1:80+80*8]
colnames(seqtab.nochim_df_check)[1+1:80+80*9]
colnames(seqtab.nochim_df_check)[1+1:80+80*10]
colnames(seqtab.nochim_df_check)[1+1:80+80*11]
colnames(seqtab.nochim_df_check)[1+1:80+80*12]

##############################################
# Summary tables
seqtab.nochim4 <- t(seqtab.nochim_df3)
asv_seqs <- colnames(seqtab.nochim4)
asv_headers <- vector(dim(seqtab.nochim4)[2], mode="character")

for (i in 1:dim(seqtab.nochim4)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, paste0("Raw-data/sequence/All_agmicrobiome/",name,"/seqtab/4seqtab_16S.fa"))

# count table:
asv_tab <- t(seqtab.nochim4)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, paste0("Raw-data/sequence/All_agmicrobiome/",name,"/seqtab/4counts.tsv"), sep="\t", quote=F, col.names=NA)


# save all objects in the environment
save(list=ls(), file=paste0("Raw-data/sequence/All_agmicrobiome/",name,"list.RData"))
# clear all objects in the environment
rm(list = ls())








###################
# Fungi
###################

### need to replace these samples with redos:
to_replace <- c("110-ITS_2021-07", 
                "118-ITS_2021-07", 
                "202-ITS_2021-07",
                "215-ITS_2021-07",
                "216-ITS_2021-07",
                "318-ITS_2021-07",
                "219-ITS_2021-09",               
                "220-ITS_2021-09")
replacement <- c("110r-ITS_2021-07", 
                 "118r-ITS_2021-07", 
                 "202r-ITS_2021-07",
                 "215r-ITS_2021-07",
                 "216r-ITS_2021-07",
                 "318r-ITS_2021-07",
                 "219_ITS_2021-fall",
                 "220_ITS_2021-fall")

# sample dates
sampdate <- c("2020-09", # fall 2020, 
              "2021-04", # spring 2021, 
              "2021-07", # summer 2021
              "2021-09", # fall 2021
              "2022-spring", # spring 2022
              "2022-summer", # summer 2022
              "2022-fall 2023-spring summer", # fall 2022, spring 2023, some of summer 2023
              "2023 summer-2024 fall") # rest of summer 2023, fall 2023, spring 2024, summer 2024, fall 2024
region <- "ITS"
name <- "Fungi"



# read in data from each sampling period (post-chimera-removal)
tab1 <- t(read.csv(paste0("Raw-data/sequence/",sampdate[1],"_agmicrobiome/",name,"/seqtab/2seqtab.csv"), row.names = 1)) 
tab2 <- t(read.csv(paste0("Raw-data/sequence/",sampdate[2],"_agmicrobiome/",name,"/seqtab/2seqtab.csv"), row.names = 1)) 
tab3 <- t(read.csv(paste0("Raw-data/sequence/",sampdate[3],"_agmicrobiome/",name,"/seqtab/2seqtab.csv"), row.names = 1)) 
tab4 <- t(read.csv(paste0("Raw-data/sequence/",sampdate[4],"_agmicrobiome/",name,"/seqtab/2seqtab.csv"), row.names = 1)) 
tab5 <- t(read.csv(paste0("Raw-data/sequence/",sampdate[5],"_agmicrobiome/",name,"/seqtab/2seqtab.csv"), row.names = 1)) 
tab6 <- t(read.csv(paste0("Raw-data/sequence/",sampdate[6],"_agmicrobiome/",name,"/seqtab/2seqtab.csv"), row.names = 1)) 
tab7 <- readRDS(file=paste0("Raw-data/sequence/2022-fall 2023-spring summer_agmicrobiome/Fungi/rocky_dada_ITS/output/02_nochimera_mergeruns/rocky_dada_ITS/rocky_dada_ITS_seqtab_final_renamed.rds"))
tab8 <- readRDS(file=paste0("Raw-data/sequence/2023 summer-2024 fall_agmicrobiome/ITS/rocky_dada_ITS/output/02_nochimera_mergeruns/rocky_dada_ITS/rocky_dada_ITS_seqtab_final_renamed.rds"))

# The sequence table is a matrix with rows corresponding to (and named by) the 
# samples, and columns corresponding to (and named by) the sequence variants. 
rownames(tab1)[1]
rownames(tab2)[1]
rownames(tab3)[1]
rownames(tab4)[1]
rownames(tab5)[1]
rownames(tab6)[1]
rownames(tab7)[1]
rownames(tab8)[1]



# clean up sample names
rownames(tab1) <- gsub("X", "", rownames(tab1))
rownames(tab1) <- gsub(".1.fungi", "-ITS", rownames(tab1)) 
rownames(tab1) <- gsub("2020.09", "2020-09", rownames(tab1)) 
rownames(tab1) <- gsub(".fungi", "-ITS", rownames(tab1)) 

rownames(tab2) <- gsub("X", "", rownames(tab2))
rownames(tab2) <- gsub(".0.10.fungal", "-ITS", rownames(tab2))
rownames(tab2) <- gsub("2021.04", "2021-04", rownames(tab2))

rownames(tab3) <- gsub("X", "", rownames(tab3))
rownames(tab3) <- gsub("_ITS", "-ITS", rownames(tab3))
rownames(tab3) <- gsub("2021.07", "2021-07", rownames(tab3))

rownames(tab4) <- gsub("X", "", rownames(tab4))
rownames(tab4) <- gsub("_ITS", "-ITS", rownames(tab4))
rownames(tab4) <- gsub("2021.09", "2021-09", rownames(tab4))

rownames(tab5) <- gsub("X", "", rownames(tab5))
rownames(tab5) <- gsub("_ITS", "-ITS", rownames(tab5))
rownames(tab5) <- gsub("2022.spring", "2022-spring", rownames(tab5))

rownames(tab6) <- gsub("X", "", rownames(tab6))
rownames(tab6) <- gsub("_ITS", "-ITS", rownames(tab6))
rownames(tab6) <- gsub("2022.summer", "2022-summer", rownames(tab6))

rownames(tab7) <- gsub("-1", "_ITS", rownames(tab7))
rownames(tab7) <- gsub("Fall_2021", "2021-fall", rownames(tab7))
rownames(tab7) <- gsub("Fall_2022", "2022-fall", rownames(tab7))
rownames(tab7) <- gsub("Spring_2023", "2023-spring", rownames(tab7))
rownames(tab7) <- gsub("Summer_2023", "2023-summer", rownames(tab7))

rownames(tab8) <- gsub("-1", "_ITS", rownames(tab8))
rownames(tab8) <- gsub("Fall_2021", "2021-fall", rownames(tab8))
rownames(tab8) <- gsub("Fall_2022", "2022-fall", rownames(tab8))
rownames(tab8) <- gsub("Spring_2023", "2023-spring", rownames(tab8))
rownames(tab8) <- gsub("Summer_2023", "2023-summer", rownames(tab8))




# inspect sequence lengths:
table(nchar(colnames(tab1)))
table(nchar(colnames(tab2)))
table(nchar(colnames(tab3)))
table(nchar(colnames(tab4)))
table(nchar(colnames(tab5)))
table(nchar(colnames(tab6)))
table(nchar(colnames(tab7)))
table(nchar(colnames(tab8)))
# Make sure the lengths of the merged sequences fall within the
# expected range for your targeted amplicon.

# NOTE: Sequences that are much longer or shorter than expected may be the 
# result of non-specific priming. You can remove non-target-length sequences 
# from your sequence table: 
# eg. seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256] 
# This is analogous to “cutting a band” in-silico to get amplicons of the 
# targeted length.

# make sure the plots are in order
rownames(tab1)
rownames(tab2)
rownames(tab3)
rownames(tab4)
rownames(tab5)
rownames(tab6)
rownames(tab7)
rownames(tab8)

tab1 <-  tab1[order(rownames(tab1)),]
tab2 <-  tab2[order(rownames(tab2)),]
tab3 <-  tab3[order(rownames(tab3)),]
tab4 <-  tab4[order(rownames(tab4)),]
tab5 <-  tab5[order(rownames(tab5)),]
tab6 <-  tab6[order(rownames(tab6)),]
#tab7 <-  tab7[order(rownames(tab7)),]
#tab8 <-  tab8[order(rownames(tab8)),]



# Combine datasets from different sampling periods:
seqtab_combo <- mergeSequenceTables(tab1, tab2, tab3, tab4, tab5, tab6, as.matrix(tab7), as.matrix(tab8))


# remove bad samples
seqtab_combo_new <- seqtab_combo
# check that we are removing the right samples
dat <- data.frame(a = rownames(seqtab_combo_new[which(rownames(seqtab_combo_new) %in% to_replace),]), 
           b = to_replace, 
           c= rownames(seqtab_combo_new[which(rownames(seqtab_combo_new) %in% replacement),]), 
           d = replacement); dat
# replace
seqtab_combo_new[which(rownames(seqtab_combo_new) %in% to_replace),] <- seqtab_combo_new[which(rownames(seqtab_combo_new) %in% replacement),]
seqtab_combo_new[which(rownames(seqtab_combo_new) %in% c(to_replace[1], replacement[1])),1:2] # check
# remove old columns for replacements
seqtab_combo_new <- seqtab_combo_new[-which(rownames(seqtab_combo_new) %in% replacement),]

seqtab_combo <- seqtab_combo_new

rownames(seqtab_combo)


##############################################
# Remove sequences that are found within the PCR blank,
# and then remove sequences with less than 10 reads

seqtab.nochim_df <- as.data.frame(t(seqtab_combo))
# determine if any asvs were detected in pcr blank
colnames(seqtab.nochim_df)
unique(seqtab.nochim_df$`PCRblank-ITS_2020-09`); which(seqtab.nochim_df$`PCRblank-ITS_2020-09`>0) 
unique(seqtab.nochim_df$`Kivlin_FungiBlank-ITS_2021-09`); which(seqtab.nochim_df$`Kivlin_FungiBlank-ITS_2021-09`>0)
unique(seqtab.nochim_df$`Blank-ITS_2022-summer`); which(seqtab.nochim_df$`Blank-ITS_2022-summer`>0)
unique(seqtab.nochim_df$`Blank-ITS_2022-spring`); which(seqtab.nochim_df$`Blank-ITS_2022-spring`>0)
unique(seqtab.nochim_df$blank1_Fall_2024); which(seqtab.nochim_df$blank1_Fall_2024>0)
unique(seqtab.nochim_df$blank2_Fall_2024); which(seqtab.nochim_df$blank2_Fall_2024>0)

# there were asvs detected in the blanks, but they have <10 reads so they  
# will be removed anyway



# remove column with pcr blank
blankcol <- which(colnames(seqtab.nochim_df) %in% c("PCRblank-ITS_2020-09", "Kivlin_FungiBlank-ITS_2021-09", "Blank-ITS_2022-summer", "Blank-ITS_2022-spring", "blank1_Fall_2024", "blank2_Fall_2024"))
seqtab.nochim_df2 <- seqtab.nochim_df[,-blankcol]



#### remove asvs with less than 10 reads across samples
seqtab.nochim_df2$sum <- rowSums(seqtab.nochim_df2)
length(seqtab.nochim_df2$sum) # check against 1seqtab-dim.csv

# number of asvs with less than 10 reads
write.csv(length(which(seqtab.nochim_df2$sum<10)), file=paste0("Raw-data/sequence/All_agmicrobiome/",name,"/seqtab/3lessthan10reads.csv"))
# number of asvs with > 10 reads
write.csv(length(which(seqtab.nochim_df2$sum>=10)), file=paste0("Raw-data/sequence/All_agmicrobiome/",name,"/seqtab/3greaterthan10reads.csv"))
# proportion of asvs retained
write.csv(length(which(seqtab.nochim_df2$sum>=10))/length(seqtab.nochim_df2$sum), file=paste0("Raw-data/sequence/All_agmicrobiome/",name,"/seqtab/3proportionretained.csv"))

# remove
seqtab.nochim_df3 <- seqtab.nochim_df2[-which(seqtab.nochim_df2$sum<10),]
hist(seqtab.nochim_df3$sum); range(seqtab.nochim_df3$sum) # check
dim(seqtab.nochim_df3)


write.csv(seqtab.nochim_df3, file=paste0("Raw-data/sequence/All_agmicrobiome/",name,"/seqtab/4seqtab_ITS.csv"))

seqtab.nochim_df3_ITS <- read.csv(file=paste0("Raw-data/sequence/All_agmicrobiome/Fungi/seqtab/4seqtab_ITS.csv"))


# check sample names
seqtab.nochim_df_check <- read.csv(file=paste0("Raw-data/sequence/All_agmicrobiome/Fungi/seqtab/4seqtab_ITS.csv"))
colnames(seqtab.nochim_df_check)[1+1:80+80*0]
colnames(seqtab.nochim_df_check)[1+1:80+80*1]
colnames(seqtab.nochim_df_check)[1+1:80+80*2]
colnames(seqtab.nochim_df_check)[1+1:80+80*3]
colnames(seqtab.nochim_df_check)[1+1:80+80*4]
colnames(seqtab.nochim_df_check)[1+1:80+80*5]
colnames(seqtab.nochim_df_check)[1+1:80+80*6]
colnames(seqtab.nochim_df_check)[1+1:80+80*7]
colnames(seqtab.nochim_df_check)[1+1:80+80*8]
colnames(seqtab.nochim_df_check)[1+1:80+80*9]
colnames(seqtab.nochim_df_check)[1+1:80+80*10]
colnames(seqtab.nochim_df_check)[1+1:80+80*11]
colnames(seqtab.nochim_df_check)[1+1:80+80*12]



##############################################
# Summary tables
seqtab.nochim4 <- t(seqtab.nochim_df3)
asv_seqs <- colnames(seqtab.nochim4)
asv_headers <- vector(dim(seqtab.nochim4)[2], mode="character")

for (i in 1:dim(seqtab.nochim4)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, paste0("Raw-data/sequence/All_agmicrobiome/",name,"/seqtab/4seqtab_ITS.fa"))

# count table:
asv_tab <- t(seqtab.nochim4)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, paste0("Raw-data/sequence/All_agmicrobiome/",name,"/seqtab/4counts.tsv"), sep="\t", quote=F, col.names=NA)


# save all objects in the environment
save(list=ls(), file=paste0("Raw-data/sequence/All_agmicrobiome/",name,"list.RData"))
# clear all objects in the environment
rm(list = ls())







