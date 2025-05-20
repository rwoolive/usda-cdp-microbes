##################################################
###### 1. calculate diversity 




library(vegan)
# install.packages("remotes")
# remotes::install_github("jfq3/QsRutils")
library(QsRutils)
library(ecodist)
library(ape)
library(lme4)
library(tidyverse)
library(emmeans)
library(nlme)
library(plyr)
library(multcomp)
library(ggpubr)
library(miaViz)
library(phyloseq)


nplots <- 80 # number of plots within the field site


# field plot metadata
plotdat <- read.csv("Raw-data/WTREC-CDP_plot-data.csv")
# field plot soil property data
soildat <- read.csv("Processed-data/WTREC-CDP_processed-data-simple.csv")
# subset soildat to 0-10cm depth
soildat <- soildat[which(soildat$Depth=="0-10 cm"),] 
soildat$season.year <- paste(soildat$Season, soildat$Year, sep=" ")
soildat$season.year <- as.factor(soildat$season.year )
soildat$season.year <- factor(soildat$season.year, levels(soildat$season.year)[c(1,6,10,2,7,11,3,8,12,4,9,13,5)])






################  Fungi ################  

## there are 13 sampling points:
# fall 2020, (baseline)
# spring 2021, summer 2021, fall 2021, 
# spring 2022, summer 2022, fall 2022
# spring 2023, summer 2023, fall 2023
# spring 2024, summer 2024, fall 2024

sampdate <- c("All_agmicrobiome") 
region <- "ITS"
name <- "Fungi"
# number of sampling periods for which we have sequence data
nsamps <- 13








# asv taxonomy
taxa <- read.csv(paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/5_taxonomy-assignment-with-ft.csv"))
# # subset to functional type
# taxa <- taxa[which(taxa$Primary_lifestyle %in% c(...))]

taxa$asv_number <- as.numeric(1:dim(taxa)[1])

# asv abundances
abundances <- read.csv(paste0("Raw-data/sequence/",sampdate,"/",name,"/seqtab/5seqtab.csv"))
abundances$asv_number <- as.numeric(abundances$X.1)
# simplify sample names
dat <- data.frame(original.name = rep(NA, nplots*nsamps),
                  new.name = rep(NA, nplots*nsamps),
                  plot = rep(NA, nplots*nsamps),
                  season.year = rep(NA, nplots*nsamps))
for(i in 1:dim(dat)[1]){
  dat$original.name[i] <- colnames(abundances)[i+2]
  dat$plot[i] <- substr(colnames(abundances)[i+2],2,4)
  if(substr(colnames(abundances)[i+2],10,16)=="2020.09"){dat$season.year[i] <- "Fall_2020"}
  if(substr(colnames(abundances)[i+2],10,16)=="2021.04"){dat$season.year[i] <- "Spring_2021"}
  if(substr(colnames(abundances)[i+2],10,16)=="2021.07"){dat$season.year[i] <- "Summer_2021"}
  if(substr(colnames(abundances)[i+2],10,16)=="2021.09"){dat$season.year[i] <- "Fall_2021"}
  if(substr(colnames(abundances)[i+2],10,16)=="2022.sp"){dat$season.year[i] <- "Spring_2022"}
  if(substr(colnames(abundances)[i+2],10,16)=="2022.su"){dat$season.year[i] <- "Summer_2022"}
  if(substr(colnames(abundances)[i+2],10,16)=="2022.fa"){dat$season.year[i] <- "Fall_2022"}
  if(substr(colnames(abundances)[i+2],10,16)=="2023.sp"){dat$season.year[i] <- "Spring_2023"}
  if(substr(colnames(abundances)[i+2],10,16)=="2023.su"){dat$season.year[i] <- "Summer_2023"}
  if(substr(colnames(abundances)[i+2],10,19)=="Fall_2023"){dat$season.year[i] <- "Fall_2023"}
  if(substr(colnames(abundances)[i+2],10,20)=="Spring_2024"){dat$season.year[i] <- "Spring_2024"}
  if(substr(colnames(abundances)[i+2],10,21)=="Summer_2024"){dat$season.year[i] <- "Summer_2024"}
  if(substr(colnames(abundances)[i+2],10,19)=="Fall_2024"){dat$season.year[i] <- "Fall_2024"}
  dat$new.name[i] <- paste0(dat$plot[i],"_",dat$season.year[i])
}


for(i in 1:dim(dat)[1]){
  dat$season[i] <- strsplit(dat$season.year[i], "_")[[1]][1]
  dat$year[i] <- strsplit(dat$season.year[i], "_")[[1]][2]
}
dat$Time <- gsub(pattern = "_", replacement = " ", x = dat$season.year)
dat$season <- as.factor(dat$season)
dat$season <- factor(dat$season, levels(dat$season)[c(2,3,1)])
# ordering <- order(dat$year, dat$season, dat$plot)
# dat <- dat[ordering,]
write.csv(dat, paste0("Processed-data/sequences/sample-metadata-", region, ".csv"))

# rename column names of abundances dataframe
for(i in 1:dim(dat)[1]){
  colnames(abundances)[i+2] <- dat$new.name[which(dat$original.name==colnames(abundances)[i+2])]
}

write.csv(abundances, paste0("Raw-data/sequence/",sampdate,"/",name,"/seqtab/5seqtab_ordered.csv"))


colnames(abundances)[1:10]
abunds <- abundances[,c(3:dim(abundances)[2])]
colnames(abunds)[1:10]
#rownames(abunds) <- abundances$X.1
abunds2 <- t(abunds)

# check:
rownames(abunds2)[1:80]
rownames(abunds2)[1:80+80*1]
rownames(abunds2)[1:80+80*2]
rownames(abunds2)[1:80+80*3]
rownames(abunds2)[1:80+80*4]
rownames(abunds2)[1:80+80*5]
rownames(abunds2)[1:80+80*6]
rownames(abunds2)[1:80+80*7] 
rownames(abunds2)[1:80+80*8] 
rownames(abunds2)[1:80+80*9] 
rownames(abunds2)[1:80+80*10] 
rownames(abunds2)[1:80+80*11] 
rownames(abunds2)[1:80+80*12] 

# rarefy to this number of reads per sample
rar <- 5000





#################################################
# create stacked bar charts for relative abundance

# first, create a phyloseq object
ASV <- otu_table(abunds2[1:1040,], taxa_are_rows = F)
rownames(dat) <- dat$new.name
dat$Time <- as.factor(dat$Time)
dat$Time <- factor(dat$Time, levels(dat$Time)[c(1,6,10,2,7,11,3,8,12,4,9,13,5)])
SAM <- sample_data(dat)
TAX <- taxa[,3:11]
TAX$Phylum <- gsub(pattern="p__", replacement="", TAX$Phylum )
rownames(TAX) <- paste0("sp", row.names(taxa))
ps <- phyloseq(otu_table=ASV, tax_table=as.matrix(TAX), sam_data=SAM)
tax_table(ps) <- as.matrix(TAX)

# convert to treesummarizedexperiment object
tse <- convertFromPhyloseq(ps)

# Apply relative transform
tser <- transformAssay(tse, method = "relabundance")

# Get top phyla
tser_p <- agglomerateByRank(tser, rank="Phylum")
top_taxa <- getTop(tser_p, top = 10, assay.type = "relabundance")


# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
phylum_renamed <- lapply(rowData(tser)$Phylum, function(x){if (x %in% top_taxa) {x} else {"Other"}})
rowData(tser)$Phylum <- as.character(phylum_renamed)

# plot by phylum
p <- plotAbundance(tser, assay.type="relabundance", group="Phylum", order.row.by="abund",
                   add_legend=T, as.relative=T, col.var = "Time") 
q <- patchwork::wrap_plots(p, ncol = 1, heights = c(0.95, 0.05))
ggpubr::ggexport(q, filename = paste0("Figures/relative abundance/", region,".png"), height=1700, width=3000, res = 300)

# plot by phylum, average by time
tser_f <- agglomerateByVariable(tser, by="cols", group="Time", update.tree=T)
p <- plotAbundance(tser_f, assay.type="relabundance", group="Phylum", order.row.by="abund",
                   add_legend=T, as.relative=T) +
  theme(axis.text.x = element_text(angle = 270,  hjust=0))
ggpubr::ggexport(p, filename = paste0("Figures/relative abundance/", region,"_averaged.png"), height=1300, width=2000, res = 300)

# Apply relative transform
tser <- transformAssay(tse, method = "relabundance")

# Get top phyla
rowData(tser)$Phylum <- rowData(tser)$Primary_lifestyle
tser_p <- agglomerateByRank(tser, rank="Phylum")
top_taxa <- getTop(tser_p, top = 10, assay.type = "relabundance")


# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
phylum_renamed <- lapply(rowData(tser)$Phylum, function(x){if (x %in% top_taxa) {x} else {"Other"}})
rowData(tser)$Phylum <- as.character(phylum_renamed)

# plot by lifestyle
p <- plotAbundance(tser, assay.type="relabundance", group="Phylum", order.row.by="abund",
                   add_legend=T, as.relative=T, col.var = "Time") 
q <- patchwork::wrap_plots(p, ncol = 1, heights = c(0.95, 0.05))
ggpubr::ggexport(q, filename = paste0("Figures/relative abundance/", region,"_lifestyle.png"), height=1700, width=3000, res = 300)

# plot by lifestyle, average by time
tser_f <- agglomerateByVariable(tser, by="cols", group="Time", update.tree=T)
p <- plotAbundance(tser_f, assay.type="relabundance", group="Phylum", order.row.by="abund",
                   add_legend=T, as.relative=T) +
  theme(axis.text.x = element_text(angle = 270,  hjust=0))
ggpubr::ggexport(p, filename = paste0("Figures/relative abundance/", region,"_lifestyle_averaged.png"), height=1300, width=2000, res = 300)





# ### create curves
# maxsamp <- 100000
# 
# rownames(abunds2[1:80,])
# fall20 <- rarecurve(abunds2[1:80,], step = 1000,  cex=0.5, tidy = TRUE)
# df <- fall20 %>% dplyr::group_by(Site) %>%  dplyr::summarise(maxSample = max(Sample, na.rm=TRUE), maxSpecies = max(Species, na.rm=TRUE))
# df <- df[order(df$maxSpecies),]
# df$Site <- substr(df$Site, start=1, stop=3)
# fall20p <- ggplot(data=fall20, aes(x=Sample, y=Species, group=Site)) + lims(x=c(0,maxsamp)) +
#   geom_label(inherit.aes = F,data = as.data.frame(df), fill=NA,aes(label = Site, x=maxSample, y=maxSpecies), na.rm = TRUE, nudge_x = 1000, hjust = 0, size=2, label.size = 0) +
#   geom_line(color="blue") + geom_vline(xintercept = rar) +
#   theme_bw() + labs(x="Number of samples", y="Number of ASVs", title="Fall 2020")
# fall20p
# 
# rownames(abunds2[81:160,])
# spr21 <- rarecurve(abunds2[81:160,], step = 1000,  cex=0.5, tidy = TRUE)
# df <- spr21 %>% dplyr::group_by(Site) %>%  dplyr::summarise(maxSample = max(Sample, na.rm=TRUE), maxSpecies = max(Species, na.rm=TRUE))
# df <- df[order(df$maxSpecies),]
# df$Site <- substr(df$Site, start=1, stop=3)
# spr21p <- ggplot(data=spr21, aes(x=Sample, y=Species, group=Site)) + lims(x=c(0,maxsamp)) +
#   geom_label(inherit.aes = F,data = as.data.frame(df), fill=NA,aes(label = Site, x=maxSample, y=maxSpecies), na.rm = TRUE, nudge_x = 1000, hjust = 0, size=2, label.size = 0) +
#   geom_line(color="blue") + geom_vline(xintercept = rar) +
#   theme_bw() + labs(x="Number of samples", y="Number of ASVs", title="Spring 2021")
# spr21p
# 
# rownames(abunds2[161:240,])
# sum21 <- rarecurve(abunds2[161:240,], step = 1000,  cex=0.5, tidy = TRUE)
# df <- sum21 %>% dplyr::group_by(Site) %>%  dplyr::summarise(maxSample = max(Sample, na.rm=TRUE), maxSpecies = max(Species, na.rm=TRUE))
# df <- df[order(df$maxSpecies),]
# df$Site <- substr(df$Site, start=1, stop=3)
# sum21p <- ggplot(data=sum21, aes(x=Sample, y=Species, group=Site)) + lims(x=c(0,maxsamp)) +
#   geom_label(inherit.aes = F,data = as.data.frame(df), fill=NA,aes(label = Site, x=maxSample, y=maxSpecies), na.rm = TRUE, nudge_x = 1000, hjust = 0, size=2, label.size = 0) +
#   geom_line(color="blue") + geom_vline(xintercept = rar) +
#   theme_bw() + labs(x="Number of samples", y="Number of ASVs", title="Summer 2021")
# sum21p
# 
# rownames(abunds2[241:320,])
# fall21 <- rarecurve(abunds2[241:320,], step = 1000,  cex=0.5, tidy = TRUE)
# df <- fall21 %>% dplyr::group_by(Site) %>%  dplyr::summarise(maxSample = max(Sample, na.rm=TRUE), maxSpecies = max(Species, na.rm=TRUE))
# df <- df[order(df$maxSpecies),]
# df$Site <- substr(df$Site, start=1, stop=3)
# fall21p <- ggplot(data=fall21, aes(x=Sample, y=Species, group=Site)) + lims(x=c(0,maxsamp)) +
#   geom_label(inherit.aes = F,data = as.data.frame(df), fill=NA,aes(label = Site, x=maxSample, y=maxSpecies), na.rm = TRUE, nudge_x = 1000, hjust = 0, size=2, label.size = 0) +
#   geom_line(color="blue") + geom_vline(xintercept = rar) +
#   theme_bw() + labs(x="Number of samples", y="Number of ASVs", title="Fall 2021")
# fall21p
# 
# rownames(abunds2[321:400,])
# spr22 <- rarecurve(abunds2[321:400,], step = 1000,  cex=0.5, tidy = TRUE)
# df <- spr22 %>% dplyr::group_by(Site) %>%  dplyr::summarise(maxSample = max(Sample, na.rm=TRUE), maxSpecies = max(Species, na.rm=TRUE))
# df <- df[order(df$maxSpecies),]
# df$Site <- substr(df$Site, start=1, stop=3)
# spr22p <- ggplot(data=spr22, aes(x=Sample, y=Species, group=Site)) + lims(x=c(0,maxsamp)) +
#   geom_label(inherit.aes = F,data = as.data.frame(df), fill=NA,aes(label = Site, x=maxSample, y=maxSpecies), na.rm = TRUE, nudge_x = 1000, hjust = 0, size=2, label.size = 0) +
#   geom_line(color="blue") + geom_vline(xintercept = rar) +
#   theme_bw() + labs(x="Number of samples", y="Number of ASVs", title="Spring 2022")
# spr22p
# 
# rownames(abunds2[401:480,])
# sum22 <- rarecurve(abunds2[401:480,], step = 1000,  cex=0.5, tidy = TRUE)
# df <- sum22 %>% dplyr::group_by(Site) %>%  dplyr::summarise(maxSample = max(Sample, na.rm=TRUE), maxSpecies = max(Species, na.rm=TRUE))
# df <- df[order(df$maxSpecies),]
# df$Site <- substr(df$Site, start=1, stop=3)
# sum22p <- ggplot(data=sum22, aes(x=Sample, y=Species, group=Site)) + lims(x=c(0,maxsamp)) +
#   geom_label(inherit.aes = F,data = as.data.frame(df), fill=NA,aes(label = Site, x=maxSample, y=maxSpecies), na.rm = TRUE, nudge_x = 1000, hjust = 0, size=2, label.size = 0) +
#   geom_line(color="blue") + geom_vline(xintercept = rar) +
#   theme_bw() + labs(x="Number of samples", y="Number of ASVs", title="Summer 2022")
# sum22p
# 
# rownames(abunds2[481:560,])
# fall22 <- rarecurve(abunds2[481:560,], step = 1000,  cex=0.5, tidy = TRUE)
# df <- fall22 %>% dplyr::group_by(Site) %>%  dplyr::summarise(maxSample = max(Sample, na.rm=TRUE), maxSpecies = max(Species, na.rm=TRUE))
# df <- df[order(df$maxSpecies),]
# df$Site <- substr(df$Site, start=1, stop=3)
# fall22p <- ggplot(data=fall22, aes(x=Sample, y=Species, group=Site)) + lims(x=c(0,maxsamp)) +
#   geom_label(inherit.aes = F,data = as.data.frame(df), fill=NA,aes(label = Site, x=maxSample, y=maxSpecies), na.rm = TRUE, nudge_x = 1000, hjust = 0, size=2, label.size = 0) +
#   geom_line(color="blue") + geom_vline(xintercept = rar) +
#   theme_bw() + labs(x="Number of samples", y="Number of ASVs", title="Fall 2022")
# fall22p
# 
# rownames(abunds2[561:640,])
# spr23 <- rarecurve(abunds2[561:640,], step = 1000,  cex=0.5, tidy = TRUE)
# df <- spr23 %>% dplyr::group_by(Site) %>%  dplyr::summarise(maxSample = max(Sample, na.rm=TRUE), maxSpecies = max(Species, na.rm=TRUE))
# df <- df[order(df$maxSpecies),]
# df$Site <- substr(df$Site, start=1, stop=3)
# spr23p <- ggplot(data=spr23, aes(x=Sample, y=Species, group=Site)) + lims(x=c(0,maxsamp)) +
#   geom_label(inherit.aes = F,data = as.data.frame(df), fill=NA,aes(label = Site, x=maxSample, y=maxSpecies), na.rm = TRUE, nudge_x = 1000, hjust = 0, size=2, label.size = 0) +
#   geom_line(color="blue") + geom_vline(xintercept = rar) +
#   theme_bw() + labs(x="Number of samples", y="Number of ASVs", title="Spring 2023")
# spr23p
# 
# rownames(abunds2[641:720,])
# sum23 <- rarecurve(abunds2[641:720,], step = 1000,  cex=0.5, tidy = TRUE)
# df <- sum23 %>% dplyr::group_by(Site) %>%  dplyr::summarise(maxSample = max(Sample, na.rm=TRUE), maxSpecies = max(Species, na.rm=TRUE))
# df <- df[order(df$maxSpecies),]
# df$Site <- substr(df$Site, start=1, stop=3)
# sum23p <- ggplot(data=sum23, aes(x=Sample, y=Species, group=Site)) + lims(x=c(0,maxsamp)) +
#   geom_label(inherit.aes = F,data = as.data.frame(df), fill=NA,aes(label = Site, x=maxSample, y=maxSpecies), na.rm = TRUE, nudge_x = 1000, hjust = 0, size=2, label.size = 0) +
#   geom_line(color="blue") + geom_vline(xintercept = rar) +
#   theme_bw() + labs(x="Number of samples", y="Number of ASVs", title="Summer 2023")
# sum23p
# 
# rownames(abunds2[721:800,])
# fall23 <- rarecurve(abunds2[721:800,], step = 1000,  cex=0.5, tidy = TRUE)
# df <- fall23 %>% dplyr::group_by(Site) %>%  dplyr::summarise(maxSample = max(Sample, na.rm=TRUE), maxSpecies = max(Species, na.rm=TRUE))
# df <- df[order(df$maxSpecies),]
# df$Site <- substr(df$Site, start=1, stop=3)
# fall23p <- ggplot(data=fall23, aes(x=Sample, y=Species, group=Site)) + lims(x=c(0,maxsamp)) +
#   geom_label(inherit.aes = F,data = as.data.frame(df), fill=NA,aes(label = Site, x=maxSample, y=maxSpecies), na.rm = TRUE, nudge_x = 1000, hjust = 0, size=2, label.size = 0) +
#   geom_line(color="blue") + geom_vline(xintercept = rar) +
#   theme_bw() + labs(x="Number of samples", y="Number of ASVs", title="Fall 2023")
# fall23p
# 
# rownames(abunds2[801:880,])
# spr24 <- rarecurve(abunds2[801:880,], step = 1000,  cex=0.5, tidy = TRUE)
# df <- spr24 %>% dplyr::group_by(Site) %>%  dplyr::summarise(maxSample = max(Sample, na.rm=TRUE), maxSpecies = max(Species, na.rm=TRUE))
# df <- df[order(df$maxSpecies),]
# df$Site <- substr(df$Site, start=1, stop=3)
# spr24p <- ggplot(data=spr24, aes(x=Sample, y=Species, group=Site)) + lims(x=c(0,maxsamp)) +
#   geom_label(inherit.aes = F,data = as.data.frame(df), fill=NA,aes(label = Site, x=maxSample, y=maxSpecies), na.rm = TRUE, nudge_x = 1000, hjust = 0, size=2, label.size = 0) +
#   geom_line(color="blue") + geom_vline(xintercept = rar) +
#   theme_bw() + labs(x="Number of samples", y="Number of ASVs", title="Spring 2024")
# spr24p
# 
# rownames(abunds2[881:960,])
# sum24 <- rarecurve(abunds2[881:960,], step = 1000,  cex=0.5, tidy = TRUE)
# df <- sum24 %>% dplyr::group_by(Site) %>%  dplyr::summarise(maxSample = max(Sample, na.rm=TRUE), maxSpecies = max(Species, na.rm=TRUE))
# df <- df[order(df$maxSpecies),]
# df$Site <- substr(df$Site, start=1, stop=3)
# sum24p <- ggplot(data=sum24, aes(x=Sample, y=Species, group=Site)) + lims(x=c(0,maxsamp)) +
#   geom_label(inherit.aes = F,data = as.data.frame(df), fill=NA,aes(label = Site, x=maxSample, y=maxSpecies), na.rm = TRUE, nudge_x = 1000, hjust = 0, size=2, label.size = 0) +
#   geom_line(color="blue") + geom_vline(xintercept = rar) +
#   theme_bw() + labs(x="Number of samples", y="Number of ASVs", title="Summer 2024")
# sum24p
# 
# rownames(abunds2[961:1040,])
# fall24 <- rarecurve(abunds2[961:1040,], step = 1000,  cex=0.5, tidy = TRUE)
# df <- fall24 %>% dplyr::group_by(Site) %>%  dplyr::summarise(maxSample = max(Sample, na.rm=TRUE), maxSpecies = max(Species, na.rm=TRUE))
# df <- df[order(df$maxSpecies),]
# df$Site <- substr(df$Site, start=1, stop=3)
# fall24p <- ggplot(data=fall24, aes(x=Sample, y=Species, group=Site)) + lims(x=c(0,maxsamp)) +
#   geom_label(inherit.aes = F,data = as.data.frame(df), fill=NA,aes(label = Site, x=maxSample, y=maxSpecies), na.rm = TRUE, nudge_x = 1000, hjust = 0, size=2, label.size = 0) +
#   geom_line(color="blue") + geom_vline(xintercept = rar) +
#   theme_bw() + labs(x="Number of samples", y="Number of ASVs", title="Fall 2024")
# fall24p
# 
# gvoid <- ggplot() + theme_void()
# ggexport(ggarrange(gvoid, gvoid, fall20p,
#                    spr21p, sum21p, fall21p,
#                    spr22p, sum22p, fall22p,
#                    spr23p, sum23p, fall23p,
#                    spr24p, sum24p, fall24p,
#                    ncol=3, nrow=5), height=6500, width=4500, res=300,
#          filename = paste0("Figures/sequences/", name, "_asv abundance curve.png"))
# 


### Because of insufficient sequence coverage, bacterial data from x samples, archaeal data from x samples, and fungal data from x samples were discarded. One sample (“X”) was removed from downstream analyses due to a large disparity in diversity and community composition from all other sequences.
#### rarefy to rar reads
sampsize <- rowSums(abunds2[1:1040,]) # gives the number of sequences for each plot
sort(sampsize)[1:20]
sort(sampsize, decreasing = T)[1:20]
remove <- c(sampsize[which(sampsize<rar)], sampsize[names(sampsize) %in% c("316_Summer_2021")]) 
write.csv(data.frame(sample=names(remove), seqs=remove), paste0("Raw-data/sequence/",sampdate,"/",name,"/seqtab/6removed.csv"))
abunds3 <- abunds2[1:1040,]
abunds3 <- abunds2[which(sampsize>rar),]
abunds3 <- abunds3[-which(rownames(abunds3) %in% c("316_Summer_2021")),]

sampsize <- rowSums(abunds3) # gives the number of sequences for each plot
sort(sampsize)[1:20] # check that sample sizes are greater than rar

abunds4 <- rrarefy(abunds3, rar)

abunds5 <- as.data.frame(t(abunds4))
abunds5$asv_number <- c(1:dim(abunds5)[1])
min(rowSums(abunds5[,1:length(sampsize)])) # some ASVs now have <10 reads
# remove ASVs that have less than 10 reads
abunds5 <- abunds5[-which(rowSums(abunds5[,1:length(sampsize)])<10),]


# merge taxa and abundance data into one df
taxabund <- merge(taxa, abunds5, 
                  by.x = "asv_number", by.y = "asv_number", 
                  all.x = FALSE, all.y = TRUE) 
rownames(taxabund) <- taxabund$asv_number

########## Relative abundance of functional groups:
amf_relabund_dat <- taxabund %>% dplyr::filter(Primary_lifestyle %in% c("arbuscular_mycorrhizal")) 
amf_relabund <- colSums(amf_relabund_dat[,c(13:1033)])/rar
plantpath_relabund_dat <- taxabund %>% dplyr::filter(Primary_lifestyle %in% c("plant_pathogen")) 
plantpath_relabund <- colSums(plantpath_relabund_dat[,c(13:1033)])/rar
saps <- c("soil_saprotroph", "wood_saprotroph", "unspecified_saprotroph", "litter_saprotroph", "dung_saprotroph")
sap_relabund_dat <- taxabund %>% dplyr::filter(Primary_lifestyle %in% saps | Secondary_lifestyle %in% saps) 
sap_relabund <- colSums(sap_relabund_dat[,c(13:1033)])/rar


write.csv(taxabund, paste0("Raw-data/sequence/",sampdate,"/",name,"/seqtab/6seqtab.csv"))

taxabund <- read.csv(paste0("Raw-data/sequence/",sampdate,"/",name,"/seqtab/6seqtab.csv"))


# export taxonomy info

# number of asvs
write.csv(dim(taxabund)[1], paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/3_number-taxa_rarefied.csv"))




# look at kingdom represented
kingdom <- as.data.frame(table(taxabund$Kingdom))
kingdom$percentage <- round((kingdom$Freq/dim(taxabund)[1])*100,2)
write.csv(kingdom,  paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/4_kingdom_rarefied.csv"))

# look at phyla represented
phyla <- as.data.frame(table(taxabund$Phylum))
phyla$percentage <- round((phyla$Freq/dim(taxabund)[1])*100,2)
write.csv(phyla,  paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/4a_phyla_rarefied.csv"))

# look at classes represented
classes <- as.data.frame(table(taxabund$Class))
classes$percentage <- round((classes$Freq/dim(taxabund)[1])*100,2)
write.csv(classes, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/4b_classes_rarefied.csv"))

# look at orders represented
orders <- as.data.frame(table(taxabund$Order))
orders$percentage <- round((orders$Freq/dim(taxabund)[1])*100,2)
write.csv(orders, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/4c_orders_rarefied.csv"))

# look at families represented
family <- as.data.frame(table(taxabund$Family))
family$percentage <- round((family$Freq/dim(taxabund)[1])*100,2)
write.csv(family, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/4d_families_rarefied.csv"))

# look at genera represented
genus <- as.data.frame(table(taxabund$Genus))
genus$percentage <- round((genus$Freq/dim(taxabund)[1])*100,2)
write.csv(genus, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/4e_genera_rarefied.csv"))

# look at genera represented
Primary_lifestyle <- as.data.frame(table(taxabund$Primary_lifestyle))
Primary_lifestyle$percentage <- round((Primary_lifestyle$Freq/dim(taxabund)[1])*100,2)
write.csv(Primary_lifestyle, paste0("Raw-data/sequence/",sampdate,"/",name,"/taxonomy/4f_Primary_lifestyle_rarefied.csv"))

# check ASV numbers match up in taxabund


# divmat: matrix with ASVs as columns and plots as rows (read counts)
plots <- which(substr(colnames(taxabund), 4, 4) == "_")[-c(1)] # columns that contain abundances
divmat <- taxabund[,plots] # matrix containing abundances with samples as columns
divmat2 <- t(divmat) # matrix containing abundances with samples as rows
comp <- cbind(colnames(divmat2),
              taxabund$asv_number)
colnames(divmat2) <- taxabund$asv_number
comp <- cbind(rownames(divmat2),
              colnames(taxabund)[plots])
rownames(divmat2) <- colnames(taxabund)[plots]
divmat <- as.data.frame(divmat2)
range(colSums(divmat)) # check that all asvs have at least 10 reads
which(colSums(divmat)<10)

# dataframe for diversity estimates and treatments
divdat <- dat
divdat$Year <- as.factor(soildat$Year)
divdat$Season <- as.factor(soildat$Season)
divdat$Season.Year <- as.factor(soildat$season.year)
divdat$Cropping.system <- as.factor(soildat$Cropping.system)
divdat$Cover <- as.factor(soildat$Cover)
divdat$Replicate <- as.factor(soildat$Rep)
divdat$Cover <- factor(divdat$Cover, levels(divdat$Cover)[c(2,1,4,5,3)])
divdat$Cropping.system <- factor(divdat$Cropping.system, levels(divdat$Cropping.system)[c(1,4,3,2)])

write.csv(divdat, paste0("Raw-data/sequence/",sampdate,"/metadata.csv"))

# divmat_pa: matrix with ASVs as columns and plots as rows (proportional abundance)
# Transform to get proportional abundance (divide by total number of reads within each sample)
divmat_pa <- proportions(as.matrix(divmat),1)
range(rowSums(divmat_pa)) # check to see that all rows sum up to 1
divmat_pa <- data.frame(divmat_pa)
divmat[1:10,1:10] 
divmat_pa[1:10,1:10] 
write.csv(divmat_pa,  paste0("Raw-data/sequence/",sampdate,"/divmat_pa_", region, ".csv"))

# # divmat_presabs: matrix with ASVs as columns and plots as rows (presence/absence)
 divmat_presabs <- divmat_pa
 divmat_presabs2 <- ifelse(divmat_presabs > 0, 1, 0)
range(rowSums(divmat_presabs2)) # check to see that all rows sum up to integers
divmat[1:10,1:10] 
divmat_presabs2[1:10,1:10] # compare to original divmat matrix


# look at summary of reads for each sample
summary(rowSums(divmat))
summary(rowSums(divmat_pa))
summary(rowSums(divmat_presabs2))


# ###################################################
# ###### 1. Determine if there is a difference in diversity by treatment.
# ### Calculate alpha diversity metrics in the Vegan package in R (Oksanen et al., 2008)
# # A. richness: number of species (no incorporation of rarity or dominance)
# # B. Shannon: diversity index based a natural logarithm, accounting for both abundance
# # and evenness of the taxa present; the higher the number, the greater the richness and
# # evenness of the species in the sample
# # C. Simpson: diversity index ranging from 0 to 1 which gives more weight to the
# # dominant species; the closer to 1 the greater the diversity
# # D. Inverse Simpson: diversity index ranging from 0 to 1 which gives more weight to the
# # dominant species; the closer to 1 the greater the diversity
# ###  Then, use ANOVA to determine differences in richness and diversity by treatments.
# ### It is a general understanding that a diversified cropping system increases
# ### microbial diversity – here we will experimentally test if that understanding
# ### holds true in the southeast U.S. cropping systems.
# 
# 
########## richness
### (total number of species)
richness <- rowSums(divmat_presabs2)
divdat$richness <- rep(NA,dim(divdat)[1])
for(i in 1:dim(divdat)[1]){
  if(divdat$new.name[i] %in% names(richness)){
    divdat$richness[i] <- richness[which(names(richness)==divdat$new.name[i])]
  }
}
hist(divdat$richness)

########## Shannon:
### (-sum p_i log(b) p_i, where p_i is the proportional abundance of species i and b is the base of the logarithm)
shannon.div <- vegan::diversity(divmat, "shannon")
divdat$shannon.div <- rep(NA,dim(divdat)[1])
for(i in 1:dim(divdat)[1]){
  if(divdat$new.name[i] %in% names(shannon.div)){
    divdat$shannon.div[i] <- shannon.div[which(names(shannon.div)==divdat$new.name[i])]
  }
}
hist(divdat$shannon.div)

########## Simpson's:
### (1 - sum p_i^2)
simpson.div <- vegan::diversity(divmat, "simpson")
divdat$simpson.div <- rep(NA,dim(divdat)[1])
for(i in 1:dim(divdat)[1]){
  if(divdat$new.name[i] %in% names(simpson.div)){
    divdat$simpson.div[i] <- simpson.div[which(names(simpson.div)==divdat$new.name[i])]
  }
}
hist(divdat$simpson.div)
divdat[which(divdat$simpson.div<0.98),] # some plots have a very low value
#divdat$simpson.div[which(divdat$simpson.div<0.99)] <- NA # remove the outlier

########## Inverse Simpson's:
### (1/sum p_i^2)
invsimpson.div <- vegan::diversity(divmat, "invsimpson")
divdat$invsimpson.div <- rep(NA,dim(divdat)[1])
for(i in 1:dim(divdat)[1]){
  if(divdat$new.name[i] %in% names(invsimpson.div)){
    divdat$invsimpson.div[i] <- invsimpson.div[which(names(invsimpson.div)==divdat$new.name[i])]
  }
}
hist(divdat$invsimpson.div)

########## Evenness:
# ***** must give this function count data, not proportional abundance
### (H/ln(total number of species))
evenness <- microbiome::evenness(t(divmat), "pielou")
divdat$evenness <- rep(NA,dim(divdat)[1])
for(i in 1:dim(divdat)[1]){
  if(divdat$new.name[i] %in% rownames(evenness)){
    divdat$evenness[i] <- evenness$pielou[which(rownames(evenness)==divdat$new.name[i])]
  }
}
hist(divdat$evenness)



######### add relative abundances

divdat <- merge(divdat, as.data.frame(amf_relabund), by.x="new.name", by.y=0, all=TRUE)
divdat <- merge(divdat, as.data.frame(plantpath_relabund), by.x="new.name", by.y=0, all=TRUE)
divdat <- merge(divdat, as.data.frame(sap_relabund), by.x="new.name", by.y=0, all=TRUE)

divdat <- divdat[order(divdat$Season.Year),]

write.csv(divdat, paste0("Processed-data/sequences/", region, "-diversity.csv"))





