

##############################################
# DADA2 is a pipeline that processes raw 16S and ITS sequence data. 
# It works entirely in R.
# It takes a fundamentally different approach than other pipelines (QIIME, mothur, UPARSE) 
# in that it doesn’t bin OTUs by 97% similarity. It doesn’t bin them at all, and instead 
# identifies “sequence variants” up to 1 nucleotide different between sequences.

# the main tutorial is here: https://benjjneb.github.io/dada2/index.html
# and the tutorial curated for analysis UT genomics core data is here: file:///Users/rachelwooliver/Library/Containers/com.apple.mail/Data/Library/Mail%20Downloads/B273FB41-2BC6-4285-8301-2C3265275660/DADA2Script_2020_12_14.html

# downloading sequence data from aws: 
# in terminal, run the following code
# aws s3 cp s3://agmicrobiome . --recursive
# the data will show up in the Users/rachelwooliver directory


sampdate <- "2022-spring"
region <- "ITS"
name <- "Fungi"
tlens <- c(250,210) # set truncation lengths (based on qscore plots), 
### NOTE quality scores are not great for this set
forprim <- "AACTTTYRRCAAYGGATCWCT" # forward primer
revprim <- "AGCCTCCGCTTATTGATATGCTTAART" # reverse primer


##############################################
# Starting point
# Curating Data from the Sequencing Core
# The UTK CEB core returns files in a separate folder for each sample with a 
# forward read (R1_001.fastq.gz) and a reverse read (R2_001.fastq.gz) for each 
# sample. These sequences are compressed .fastq files and can stay in this 
# format for DADA2. Fastq files are DNA sequences with associated quality 
# scores for each base in that sequence.
# First, place filtered files in filtered/ subdirectory




##############################################
# install DADA2 from Bioconductor
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("dada2", version )
# Bioconductor version 3.12 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)
# or
# install from source
# install.packages("devtools")
# library("devtools")
# devtools::install_github("benjjneb/dada2", ref="v1.16") # change the ref argument to get other versions



##############################################
# load libraries
library(dada2)
packageVersion("dada2") # [1] ‘1.16.0’



##############################################
# Directory containing the fastq files after unzipping.
path <- paste0("Raw-data/sequence/",sampdate,"_agmicrobiome/", name)
list.files(path)

##############################################
#write.csv(list.files(paste0(path,"/raw")), paste0("Raw-data/sequence/",sampdate,"_agmicrobiome/file-names0.csv"))
# Table to translate sample names
samps <- read.csv(paste0("Raw-data/sequence/",sampdate,"_agmicrobiome/file-names.csv"))
samps <- samps[which(samps$Primer==name),]


##############################################
# Read in the names of the fastq files, and perform some string manipulation 
# to get matched lists of the forward and reverse fastq files.

fnFs <- sort(list.files(paste0(path,"/raw"), pattern="_R1_001.fastq", full.names = TRUE)) # forward reads
fnRs <- sort(list.files(paste0(path,"/raw"), pattern="_R2_001.fastq", full.names = TRUE)) # reverse reads
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names0 <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#sample.names0 <- gsub("-", "_", sample.names0)

sample.names0b <- rep(NA,length(sample.names0))
for(i in 1:length(sample.names0b)){
  sample.names0b[i] <- samps$samplenames[which(samps$X.==sample.names0[i])]
}
sample.names0b <-  sapply(strsplit(sample.names0b, "-"), `[`, 1)
#sample.names0b[81] <- sample.names0[81]


sample.names <- paste0(sample.names0b,"_",region, "_",sampdate)

sample.names[1:5]
#sample.names[81]



##############################################
# Inspect read quality profiles
# We start by visualizing the quality profiles of the forward reads:

pdf(paste0("Figures/sequences/",sampdate,"/qscores_",region,"/F_1-10.pdf"), height=8, width=12)
dada2::plotQualityProfile(fnFs[c(1:10)]) # forward reads
dev.off()

pdf(paste0("Figures/sequences/",sampdate,"/qscores_",region,"/F_11-20.pdf"), height=8, width=12)
dada2::plotQualityProfile(fnFs[c(11:20)]) # forward reads
dev.off()

pdf(paste0("Figures/sequences/",sampdate,"/qscores_",region,"/F_21-30.pdf"), height=8, width=12)
dada2::plotQualityProfile(fnFs[c(21:30)]) # forward reads
dev.off()

pdf(paste0("Figures/sequences/",sampdate,"/qscores_",region,"/F_31-40.pdf"), height=8, width=12)
dada2::plotQualityProfile(fnFs[c(31:40)]) # forward reads
dev.off()

pdf(paste0("Figures/sequences/",sampdate,"/qscores_",region,"/F_41-50.pdf"), height=8, width=12)
dada2::plotQualityProfile(fnFs[c(41:50)]) # forward reads
dev.off()

pdf(paste0("Figures/sequences/",sampdate,"/qscores_",region,"/F_51-60.pdf"), height=8, width=12)
dada2::plotQualityProfile(fnFs[c(51:60)]) # forward reads
dev.off()

pdf(paste0("Figures/sequences/",sampdate,"/qscores_",region,"/F_61-70.pdf"), height=8, width=12)
dada2::plotQualityProfile(fnFs[c(61:70)]) # forward reads
dev.off()

pdf(paste0("Figures/sequences/",sampdate,"/qscores_",region,"/F_71-81.pdf"), height=8, width=12)
dada2::plotQualityProfile(fnFs[c(71:81)]) # forward reads
dev.off()



# In gray-scale is a heat map of the frequency of each quality score at 
# each base position. The mean quality score at each position is shown by 
# the green line, and the quartiles of the quality score distribution by 
# the orange lines. The red line shows the scaled proportion of reads that 
# extend to at least that position (this is more useful for other sequencing 
# technologies, as Illumina reads are typically all the same length, hence 
# the flat red line).

# DADA2 generally advises trimming the last few nucleotides to avoid less 
# well-controlled errors that can arise there. These quality profiles do 
# not suggest that any additional trimming is needed. We are looking for
# the point at which the average quality score (teal line) drops below 30.



# Now visualize the quality profile of the reverse reads:


pdf(paste0("Figures/sequences/",sampdate,"/qscores_",region,"/R_1-10.pdf"), height=8, width=12)
dada2::plotQualityProfile(fnRs[c(1:10)]) # forward reads
dev.off()

pdf(paste0("Figures/sequences/",sampdate,"/qscores_",region,"/R_11-20.pdf"), height=8, width=12)
dada2::plotQualityProfile(fnRs[c(11:20)]) # forward reads
dev.off()

pdf(paste0("Figures/sequences/",sampdate,"/qscores_",region,"/R_21-30.pdf"), height=8, width=12)
dada2::plotQualityProfile(fnRs[c(21:30)]) # forward reads
dev.off()

pdf(paste0("Figures/sequences/",sampdate,"/qscores_",region,"/R_31-40.pdf"), height=8, width=12)
dada2::plotQualityProfile(fnRs[c(31:40)]) # forward reads
dev.off()

pdf(paste0("Figures/sequences/",sampdate,"/qscores_",region,"/R_41-50.pdf"), height=8, width=12)
dada2::plotQualityProfile(fnRs[c(41:50)]) # forward reads
dev.off()

pdf(paste0("Figures/sequences/",sampdate,"/qscores_",region,"/R_51-60.pdf"), height=8, width=12)
dada2::plotQualityProfile(fnRs[c(51:60)]) # forward reads
dev.off()

pdf(paste0("Figures/sequences/",sampdate,"/qscores_",region,"/R_61-70.pdf"), height=8, width=12)
dada2::plotQualityProfile(fnRs[c(61:70)]) # forward reads
dev.off()

pdf(paste0("Figures/sequences/",sampdate,"/qscores_",region,"/R_71-81.pdf"), height=8, width=12)
dada2::plotQualityProfile(fnRs[c(71:81)]) # forward reads
dev.off()


# The reverse reads are usually of significantly worse quality, especially at the 
# end, which is common in Illumina sequencing. This isn’t too worrisome, as DADA2 
# incorporates quality information into its error model which makes the algorithm 
# robust to lower quality sequence, but trimming as the average qualities crash 
# will improve the algorithm’s sensitivity to rare sequence variants. Tell 
# DADA2 to truncate the reverse reads at a position where the quality
# distribution crashes. 
# Reverse reads always are poorer quality than forward reads. This may be optimized 
# with MiSeq v2 or v3 chemistry. My runs from the UTK CEB use v3 (technically 2x300b 
# runs) truncated at 250b to improve quality over v2 (2x250b runs). Again, each of 
# these plots is one sample. We are looking for the point at which the average 
# quality score (teal line) drops below 30.




##############################################
# Filter and trim
# Assign the filenames for the filtered fastq.gz files.

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# The goal here is to set the parameters for sequence inclusion in the final dataset. 
# We will 
# (1)Set the maximum number of expected errors per read (maxEE). The current script 
# is for 2 errors for forward reads and 5 for reverse reads. 
# (2) Remove phiX standards 
# (3) Truncate forward and reverse reads based on quality scores (truncLen). The 
# script currently is set to do this at 250b for forward reads (no truncation) and 
# 210 bases for reverse reads. 
# (4) Set the maximum number of Ns allows (default 0) The output will tell you how 
# many reads were in the original file (reads.in) and how many now remain after 
# filtering (reads.out).
# (5) Trim primers using the trimLeft argument. 


# truncQ=2 (Truncate reads at the first instance of a quality score less than or equal to truncQ)
# maxEE: 2 errors for forward reads and 5 for reverse reads


out <- dada2::filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=tlens, 
                            trimLeft = c(nchar(forprim), nchar(revprim)),
                            maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                            compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

# # If getting error when running filterAnd Trim: https://github.com/benjjneb/dada2/issues/212#issuecomment-836757295
# remove.packages("Matrix")
# # clear workspace, then restart session
# devtools::install_version("Matrix", version = "1.3.2", repos = "http://cran.us.r-project.org")
# packageVersion("Matrix") # should be 1.3.2


  
# NOTE: If you are using a less-overlapping primer set compared to V4, like V1-V2 or 
# V3-V4, your truncLen must be large enough to maintain 20 + biological.length.variation 
# nucleotides of overlap between them.

# NOTE: The standard filtering parameters are starting points, not set in stone. 
# If you want to speed up downstream computation, consider tightening maxEE. If 
# too few reads are passing the filter, consider relaxing maxEE, perhaps especially 
# on the reverse reads (eg. maxEE=c(2,5)), and reducing the truncLen to remove 
# low quality tails. Remember though, when choosing truncLen for paired-end reads 
# you must maintain overlap after truncation in order to merge them later.

# NOTE: For ITS sequencing, it is usually undesirable to truncate reads to a fixed 
# length due to the large length variation at that locus. That is OK, just leave out 
# truncLen. See the DADA2 ITS workflow for more information




##############################################
# Learn the Error Rates
# The DADA2 algorithm makes use of a parametric error model (err) and every amplicon 
# dataset has a different set of error rates. The learnErrors method learns this 
# error model from the data, by alternating estimation of the error rates and 
# inference of sample composition until they converge on a jointly consistent 
# solution. As in many machine-learning problems, the algorithm must begin with an 
# initial guess, for which the maximum possible error rates in this data are used 
# (the error rates if only the most abundant sequence is correct and all the rest 
# are errors).
# DADA2 is defined by this step. The algorithm will learn the sequence error rates 
# and correct for these going forward. This is an important step since sequences will 
# only group together if they are exactly identical. Error rates vary among datasets 
# due to sequencing region, PCR conditions, and the sequencer.



# The following runs in about 3 minutes on a 2013 Macbook Pro:
  
errF <- dada2::learnErrors(filtFs, multithread=TRUE)

errR <- dada2::learnErrors(filtRs, multithread=TRUE)

# It is always worthwhile, as a sanity check if nothing else, to visualize the 
# estimated error rates:
  
errorplot <- dada2::plotErrors(errF, nominalQ=TRUE)
pdf(paste0("Figures/sequences/",sampdate,"/errorplot_",region,"/F.pdf"), height=8, width=8)
plot(errorplot)
dev.off()

errorplot <- dada2::plotErrors(errR, nominalQ=TRUE)
pdf(paste0("Figures/sequences/",sampdate,"/errorplot_",region,"/R.pdf"), height=8, width=8)
plot(errorplot)
dev.off()


# We are looking for the points to align with the black line which is the error 
# rate that will be applied in the dada algorithm.

# The error rates for each possible transition (A→C, A→G, …) are shown. Points are 
# the observed error rates for each consensus quality score. The black line shows 
# the estimated error rates after convergence of the machine-learning algorithm. 
# The red line shows the error rates expected under the nominal definition of the 
# Q-score. Here the estimated error rates (black line) are a good fit to the 
# observed rates (points), and the error rates drop with increased quality as 
# expected. Everything looks reasonable and we proceed with confidence.




##############################################
# Dereplicate
# This speeds up the computation by merging identical sequences.


derep_forward <- derepFastq(filtFs, verbose=TRUE)
names(derep_forward) <- sample.names

derep_reverse <- derepFastq(filtRs, verbose=TRUE)
names(derep_reverse) <- sample.names



##############################################
# Sample Inference
# We are now ready to apply the core sample inference algorithm to the filtered and 
# trimmed sequence data.
# This algorithm will tell us how many unique sequences are in each sample after 
# controlling for sequencing errors. If we are interested in low abundance variants, 
# we can pool across all samples with pool=TRUE (computationally consuming).


dadaFs <- dada2::dada(derep_forward, err=errF, multithread=TRUE)

dadaRs <- dada2::dada(derep_reverse, err=errR, multithread=TRUE)

# Inspecting the returned dada-class object:
  
dadaFs[[1]]

# NOTE: DADA2 also supports 454 and Ion Torrent data, but we recommend some minor 
# parameter changes for those pyrosequencing technologies. The adventurous can 
# explore ?setDadaOpt for other adjustable algorithm parameters.

# NOTE: Extensions: By default, the dada function processes each sample 
# independently. However, pooling information across samples can increase 
# sensitivity to sequence variants that may be present at very low frequencies 
# in multiple samples. The dada2 package offers two types of pooling. 
# dada(..., pool=TRUE) performs standard pooled processing, in which all 
# samples are pooled together for sample inference. dada(..., pool="pseudo") 
# performs pseudo-pooling, in which samples are processed independently after 
# sharing information between samples, approximating pooled sample inference 
# in linear time.




##############################################
# Merge paired reads
# We now merge the forward and reverse reads together to obtain the full denoised 
# sequences. Merging is performed by aligning the denoised forward reads with the 
# reverse-complement of the corresponding denoised reverse reads, and then 
# constructing the merged “contig” sequences. By default, merged sequences are 
# only output if the forward and reverse reads overlap by at least 12 bases, and 
# are identical to each other in the overlap region (but these conditions can be 
# changed via function arguments).

# We can run this multiple ways with different amounts of overlap amount forward 
# and reverse reads. Most sequence regions for bacterial reads will overlap 
# whereas most for fungi will not. When forward and reverse reads do not overlap, 
# we can still use the data by concatenating the reads with the justConcatenate 
# function.



mergers <- dada2::mergePairs(dadaFs, derep_forward, 
                             dadaRs, derep_reverse, 
                             verbose=FALSE, justConcatenate = TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])
sink(paste0("Raw-data/sequence/",sampdate,"_agmicrobiome/",name,"/merged-reads.csv"))
print(mergers)
sink()

# Can blast some of the sequences to double check they are what you think they are


# The mergers object is a list of data.frames from each sample. Each data.frame 
# contains the merged $sequence, its $abundance, and the indices of the $forward 
# and $reverse sequence variants that were merged. Paired reads that did not 
# exactly overlap were removed by mergePairs, further reducing spurious output.

# NOTE: Considerations for your own data: Most of your reads should successfully 
# merge. If that is not the case upstream parameters may need to be revisited: 
# Did you trim away the overlap between your reads?

# NOTE: Extensions: Non-overlapping reads are supported, but not recommended, 
# with mergePairs(..., justConcatenate=TRUE).




##############################################
# Construct sequence table
# We can now construct an amplicon sequence variant table (ASV) table, a 
# higher-resolution version of the OTU table produced by traditional methods.

seqtab <- dada2::makeSequenceTable(mergers)
# row number is your sample number, column number is your ASV number
write.csv(dim(seqtab), paste0("Raw-data/sequence/",sampdate,"_agmicrobiome/",name,"/seqtab/0seqtab-dim.csv"))



# Inspect distribution of sequence lengths
write.csv(table(nchar(dada2::getSequences(seqtab))), paste0("Raw-data/sequence/",sampdate,"_agmicrobiome/",name,"/seqtab/0seqtab-sequence-length-distribution.csv"))
# Make sure the lengths of the merged sequences fall within the
# expected range for your targeted amplicon.


# If this run is part of a big data project, save your seqtab as a rds file
saveRDS(seqtab,paste0("Raw-data/sequence/",sampdate,"_agmicrobiome/",name,"/seqtab/0seqtab.csv")) 

# The sequence table is a matrix with rows corresponding to (and named by) the 
# samples, and columns corresponding to (and named by) the sequence variants. 


# NOTE: Sequences that are much longer or shorter than expected may be the 
# result of non-specific priming. You can remove non-target-length sequences 
# from your sequence table: 
# eg. seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256] 
# This is analogous to “cutting a band” in-silico to get amplicons of the 
# targeted length.



##############################################
# Remove chimeras
# The core dada method corrects substitution and indel errors, but chimeras remain. 
# Fortunately, the accuracy of sequence variants after denoising makes identifying 
# chimeric ASVs simpler than when dealing with fuzzy OTUs. Chimeric sequences are 
# identified if they can be exactly reconstructed by combining a left-segment and 
# a right-segment from two more abundant “parent” sequences.

seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

write.csv(dim(seqtab.nochim), paste0("Raw-data/sequence/",sampdate,"_agmicrobiome/",name,"/seqtab/1seqtab-dim.csv"))


chimtab <- data.frame(name = c("asvs.that.are.chimeric", "sequences.that.are.chimeric"),
                      proportion = c(1-dim(seqtab.nochim)[2]/dim(seqtab)[2], 
                                     1-sum(seqtab.nochim)/sum(seqtab)))
write.csv(chimtab, paste0("Raw-data/sequence/",sampdate,"_agmicrobiome/",name,"/seqtab/1seqtab-chimeras.csv"))


# The frequency of chimeric sequences varies substantially from dataset to dataset, 
# and depends on on factors including experimental procedures and sample complexity. 
# Here chimeras make up about X% of the merged sequence variants, but when we 
# account for the abundances of those variants we see they account for only about 
# X% of the merged sequence reads.

# NOTE: Most of your reads should remain after chimera removal (it is not uncommon 
# for a majority of sequence variants to be removed though). If most of your reads 
# were removed as chimeric, upstream processing may need to be revisited. In almost 
# all cases this is caused by primer sequences with ambiguous nucleotides that were 
# not removed prior to beginning the DADA2 pipeline.





##############################################
# Track reads through the pipeline
# As a final check of our progress, we’ll look at the number of reads that made 
# it through each step in the pipeline:
  
getN <- function(x) sum(dada2::getUniques(x))
track <- cbind(out, 
               sapply(dadaFs, getN), 
               sapply(dadaRs, getN), 
               sapply(mergers, getN), 
               rowSums(seqtab.nochim),
               round(rowSums(seqtab.nochim)/out[,1]*100, 1))
colnames(track) <- c("input", 
                     "filtered", 
                     "denoisedF", 
                     "denoisedR", 
                     "merged", 
                     "nonchim",
                     "final_perc_reads_retained")
rownames(track) <- sample.names
head(track)
track <- as.data.frame(track)
write.csv(track, paste0("Raw-data/sequence/",sampdate,"_agmicrobiome/",name,"/track.csv"))

# If processing a single sample, remove the sapply calls: e.g. replace 
# sapply(dadaFs, getN) with getN(dadaFs)


# Looks good! We kept the majority of our raw reads, and there is no over-large 
# drop associated with any single step.

# NOTE: This is a great place to do a last sanity check. Outside of filtering, 
# there should no step in which a majority of reads are lost. If a majority of 
# reads failed to merge, you may need to revisit the truncLen parameter used in 
# the filtering step and make sure that the truncated reads span your amplicon. 
# If a majority of reads were removed as chimeric, you may need to revisit the 
# removal of primers, as the ambiguous nucleotides in unremoved primers interfere 
# with chimera identification.


##############################################
# Save the files
seqtab.nochim <- t(seqtab.nochim)
write.csv(seqtab.nochim, file=paste0("Raw-data/sequence/",sampdate,"_agmicrobiome/",name,"/seqtab/2seqtab.csv"))
save(list=ls(), file=paste0(path,"list.RData"))




##############################################
# Summary tables
seqtab.nochim <- t(seqtab.nochim)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, paste0("Raw-data/sequence/",sampdate,"_agmicrobiome/",name,"/seqtab/2seqtab.fa"))

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, paste0("Raw-data/sequence/",sampdate,"_agmicrobiome/",name,"/seqtab/2counts.tsv"), sep="\t", quote=F, col.names=NA)




