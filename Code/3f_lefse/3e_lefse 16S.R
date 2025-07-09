
rm(list = ls())

#### 3e. LEfSe: identified taxonomic groups whose abundances were 
# associated with specific treatments 


# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# 
# BiocManager::install("microbiomeMarker")

library(microbiomeMarker)
library(tidyr)
library(ggplot2)
library(phyloseq)
library(ggtree)
library(stringr)
library(dplyr)
library(ape)


cropsyscols <- c('Corn-Corn-Corn' = "darkgoldenrod4", 'Soybean-Soybean-Soybean' = "darkolivegreen", 'Corn-Soybean-Corn' = "darkolivegreen3", 'Corn-Cotton-Soybean' = "deepskyblue3")
covercols <- c('No cover' = "chocolate", 'Wheat' = "goldenrod", 'Clover' = "forestgreen", 'Wheat-Clover' = "darkcyan", 'SHM' = "darkorchid4")

reso <- 400
w_cutoff <- 0.01
kw_cutoff <- 0.01
l_cutoff <- 4


#### load in all the data
# field plot metadata
plotdat <- read.csv("Raw-data/WTREC-CDP_plot-data.csv")
# field plot soil property data
soildat <- read.csv("Processed-data/WTREC-CDP_processed-data-simple.csv")
# subset soildat to 0-10cm depth
soildat <- soildat[which(soildat$Depth=="0-10 cm"),] 
soildat$season.year <- paste(soildat$Season, soildat$Year, sep=" ")
soildat$season.year <- as.factor(soildat$season.year )
soildat$season.year <- factor(soildat$season.year, levels(soildat$season.year)[c(1,6,10,2,7,11,3,8,12,4,9,13,5)])

sampdate <- c("All_agmicrobiome") 
region <- "16S"
name <- "Bacteria"
# number of sampling periods for which we have sequence data
nsamps <- 13

# import data for abundances and taxonomy
taxabund <- read.csv(paste0("Raw-data/sequence/",sampdate,"/",name,"/seqtab/6seqtab.csv"))
# import sequence metadata
seqdat <- read.csv(paste0("Processed-data/sequences/16S-diversity.csv"))

### make a phyloseq object with:
# asv_table: samples as columns, asvs as rows
asv_table <- taxabund[, 10:992]
row.names(asv_table) <- taxabund$asv_number
asv_table <- phyloseq::otu_table(asv_table, taxa_are_rows = T)
# tax: clades as columns, asvs as rows
tax <- taxabund[, 4:9]
row.names(tax) <- taxabund$asv_number
tax <- phyloseq::tax_table(as.matrix(tax))
# sam: metadata as columns, samples as rows
sam <- na.omit(seqdat)
row.names(sam) <- paste0("X", na.omit(seqdat)$new.name)
sam <- phyloseq::sample_data(sam)
# phy_tree, and refseq: NA
divexp_phy <- phyloseq::phyloseq(asv_table, tax, sam)


# ### LEfSe example
# data(kostic_crc)
# kostic_crc_small <- phyloseq::subset_taxa(
#   kostic_crc,
#   Phylum == "Firmicutes"
# )
# mm_lefse <- run_lefse(
#   kostic_crc_small,
#   wilcoxon_cutoff = 0.01,
#   group = "DIAGNOSIS",
#   kw_cutoff = 0.01,
#   multigrp_strat = TRUE,
#   lda_cutoff = 4
# )
# plot_cladogram(mm_lefse, color = c("darkgreen", "red"))





### LEfSe with all experimental timepoints included
'%notin%' <- Negate('%in%')
# subset to timepoint
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %notin% c("Fall 2020")
)
# lefse analysis
mm_lefse_0 <- run_lefse(divexp_phy_small, group = "Cover", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

# which treatments are represented?
gps <- unique(mm_lefse_0@marker_table$enrich_group)

# plot with cladogram
p_overall <- plot_cladogram(mm_lefse_0, color = covercols[gps], 
                      alpha = 0.5, clade_label_level = 5) +
  theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none")  

ggpubr::ggexport(p_overall + labs(tag = "A) Bacteria"), width = 4000, height = 2500, res = reso, 
                 filename = paste0("Figures/microbial-diversity/lefse/", region, "_overall_cover.png"))





### Corn-Corn-Corn LEfSe 
'%notin%' <- Negate('%in%')
# subset to timepoint
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %notin% c("Fall 2020") & Cropping.system %in% "Corn-Corn-Corn"
)
# lefse analysis
mm_lefse_0 <- run_lefse(divexp_phy_small, group = "Cover", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

# # which treatments are represented?
# gps <- unique(mm_lefse_0@marker_table$enrich_group)
# 
# # plot with cladogram
# p_overall <- plot_cladogram(mm_lefse_0, color = covercols[gps], 
#                             alpha = 0.5, clade_label_level = 5) +
#   theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none")  
# 
# ggpubr::ggexport(p_overall + labs(tag = "Bacteria"), width = 4000, height = 2500, res = reso, 
#                  filename = paste0("Figures/microbial-diversity/lefse/", region, "_Corn-Corn-Corn_cover.png"))

### Corn-Soybean-Corn LEfSe 
'%notin%' <- Negate('%in%')
# subset to timepoint
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %notin% c("Fall 2020") & Cropping.system %in% "Corn-Soybean-Corn"
)
# lefse analysis
mm_lefse_0 <- run_lefse(divexp_phy_small, group = "Cover", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

# # which treatments are represented?
# gps <- unique(mm_lefse_0@marker_table$enrich_group)
# 
# # plot with cladogram
# p_overall <- plot_cladogram(mm_lefse_0, color = covercols[gps], 
#                             alpha = 0.5, clade_label_level = 5) +
#   theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none")  
# 
# ggpubr::ggexport(p_overall + labs(tag = "Bacteria"), width = 4000, height = 2500, res = reso, 
#                  filename = paste0("Figures/microbial-diversity/lefse/", region, "_Corn-Soybean-Corn_cover.png"))

### Soybean-Soybean-Soybean LEfSe 
'%notin%' <- Negate('%in%')
# subset to timepoint
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %notin% c("Fall 2020") & Cropping.system %in% "Soybean-Soybean-Soybean"
)
# lefse analysis
mm_lefse_0 <- run_lefse(divexp_phy_small, group = "Cover", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

# # which treatments are represented?
# gps <- unique(mm_lefse_0@marker_table$enrich_group)
# 
# # plot with cladogram
# p_overall <- plot_cladogram(mm_lefse_0, color = covercols[gps], 
#                             alpha = 0.5, clade_label_level = 5) +
#   theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none")  
# 
# ggpubr::ggexport(p_overall + labs(tag = "Bacteria"), width = 4000, height = 2500, res = reso, 
#                  filename = paste0("Figures/microbial-diversity/lefse/", region, "_Soybean-Soybean-Soybean_cover.png"))

### Corn-Cotton-Soybean LEfSe 
'%notin%' <- Negate('%in%')
# subset to timepoint
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %notin% c("Fall 2020") & Cropping.system %in% "Corn-Cotton-Soybean"
)
# lefse analysis
mm_lefse_0 <- run_lefse(divexp_phy_small, group = "Cover", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

# # which treatments are represented?
# gps <- unique(mm_lefse_0@marker_table$enrich_group)
# 
# # plot with cladogram
# p_overall <- plot_cladogram(mm_lefse_0, color = covercols[gps], 
#                             alpha = 0.5, clade_label_level = 5) +
#   theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none")  
# 
# ggpubr::ggexport(p_overall + labs(tag = "Bacteria"), width = 4000, height = 2500, res = reso, 
#                  filename = paste0("Figures/microbial-diversity/lefse/", region, "_Corn-Cotton-Soybean_cover.png"))

### No cover LEfSe 
'%notin%' <- Negate('%in%')
# subset to timepoint
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %notin% c("Fall 2020", "Spring 2021") & Cover %in% "No cover"
)
# lefse analysis
mm_lefse_0 <- run_lefse(divexp_phy_small, group = "Cropping.system", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

# # which treatments are represented?
# gps <- unique(mm_lefse_0@marker_table$enrich_group)
# 
# # plot with cladogram
# p_overall <- plot_cladogram(mm_lefse_0, color = cropsyscols[gps],
#                             alpha = 0.5, clade_label_level = 5) +
#   theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none")
# 
# ggpubr::ggexport(p_overall + labs(tag = "Fungi"), width = 4000, height = 2500, res = reso,
#                  filename = paste0("Figures/microbial-diversity/lefse/", region, "_No cover_cropsys.png"))

### Wheat LEfSe 
'%notin%' <- Negate('%in%')
# subset to timepoint
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %notin% c("Fall 2020", "Spring 2021") & Cover %in% "Wheat"
)
# lefse analysis
mm_lefse_0 <- run_lefse(divexp_phy_small, group = "Cropping.system", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

# # which treatments are represented?
# gps <- unique(mm_lefse_0@marker_table$enrich_group)
# 
# # plot with cladogram
# p_overall <- plot_cladogram(mm_lefse_0, color = cropsyscols[gps],
#                             alpha = 0.5, clade_label_level = 5) +
#   theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none")
# 
# ggpubr::ggexport(p_overall + labs(tag = "Fungi"), width = 4000, height = 2500, res = reso,
#                  filename = paste0("Figures/microbial-diversity/lefse/", region, "_Wheat_cropsys.png"))


### Clover LEfSe 
'%notin%' <- Negate('%in%')
# subset to timepoint
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %notin% c("Fall 2020", "Spring 2021") & Cover %in% "Clover"
)
# lefse analysis
mm_lefse_0 <- run_lefse(divexp_phy_small, group = "Cropping.system", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

# # which treatments are represented?
# gps <- unique(mm_lefse_0@marker_table$enrich_group)
# 
# # plot with cladogram
# p_overall <- plot_cladogram(mm_lefse_0, color = cropsyscols[gps],
#                             alpha = 0.5, clade_label_level = 5) +
#   theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none")
# 
# ggpubr::ggexport(p_overall + labs(tag = "Fungi"), width = 4000, height = 2500, res = reso,
#                  filename = paste0("Figures/microbial-diversity/lefse/", region, "_Clover_cropsys.png"))




### Wheat-clover mix LEfSe 
'%notin%' <- Negate('%in%')
# subset to timepoint
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %notin% c("Fall 2020", "Spring 2021") & Cover %in% "Wheat-Clover"
)
# lefse analysis
mm_lefse_0 <- run_lefse(divexp_phy_small, group = "Cropping.system", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

# # which treatments are represented?
# gps <- unique(mm_lefse_0@marker_table$enrich_group)
# 
# # plot with cladogram
# p_overall <- plot_cladogram(mm_lefse_0, color = cropsyscols[gps],
#                             alpha = 0.5, clade_label_level = 5) +
#   theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none")
# 
# ggpubr::ggexport(p_overall + labs(tag = "Fungi"), width = 4000, height = 2500, res = reso,
#                  filename = paste0("Figures/microbial-diversity/lefse/", region, "_Wheat-Clover_cropsys.png"))




### SHM LEfSe 
'%notin%' <- Negate('%in%')
# subset to timepoint
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %notin% c("Fall 2020", "Spring 2021") & Cover %in% "SHM"
)
# lefse analysis
mm_lefse_0 <- run_lefse(divexp_phy_small, group = "Cropping.system", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

# # which treatments are represented?
# gps <- unique(mm_lefse_0@marker_table$enrich_group)
# 
# # plot with cladogram
# p_overall <- plot_cladogram(mm_lefse_0, color = cropsyscols[gps],
#                             alpha = 0.5, clade_label_level = 5) +
#   theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none")
# 
# ggpubr::ggexport(p_overall + labs(tag = "Fungi"), width = 4000, height = 2500, res = reso,
#                  filename = paste0("Figures/microbial-diversity/lefse/", region, "_SHM_cropsys.png"))















### Spring 2021 LEfSe
# subset to timepoint
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %in% c("Spring 2021")
)

# lefse analysis
mm_lefse_1 <- run_lefse(divexp_phy_small, group = "Cover", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

# which treatments are represented?
gps <- unique(mm_lefse_1@marker_table$enrich_group)

# plot with cladogram
p_1 <- plot_cladogram(mm_lefse_1, color = covercols[gps], 
                      alpha = 0.5, clade_label_level = 5) + 
  theme(plot.margin = margin(0, 0, 0, 0)) + 
  guides(fill = "none") #, shape = guide_legend(override.aes = list(ncol=1, color = covercols[mm_lefse_1@marker_table$enrich_group][-1], shape = 15, size = 5)))    


### Summer 2021 LEfSe
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %in% c("Summer 2021")
)
mm_lefse_2 <- run_lefse(divexp_phy_small, group = "Cover", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

gps <- unique(mm_lefse_2@marker_table$enrich_group)

p_2 <- plot_cladogram(mm_lefse_2, color = covercols[gps], 
                      alpha = 0.5, clade_label_level = 5) + 
  theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none") 


### Fall 2021 LEfSe
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %in% c("Fall 2021")
)
mm_lefse_3 <- run_lefse(divexp_phy_small, group = "Cover", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

gps <- unique(mm_lefse_3@marker_table$enrich_group)

p_3 <- plot_cladogram(mm_lefse_3, color = covercols[gps], 
                      alpha = 0.5, clade_label_level = 5) + 
  theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none") 


### Spring 2022 LEfSe
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %in% c("Spring 2022")
)
mm_lefse_4 <- run_lefse(divexp_phy_small, group = "Cover", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

gps <- unique(mm_lefse_4@marker_table$enrich_group)

p_4 <- plot_cladogram(mm_lefse_4, color = covercols[gps], 
                      alpha = 0.5, clade_label_level = 5) + 
  theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none") 


### Summer 2022 LEfSe
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %in% c("Summer 2022")
)
mm_lefse_5 <- run_lefse(divexp_phy_small, group = "Cover", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

gps <- unique(mm_lefse_5@marker_table$enrich_group)

p_5 <- plot_cladogram(mm_lefse_5, color = covercols[gps], 
                      alpha = 0.5, clade_label_level = 5) + 
  theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none") 


### Fall 2022 LEfSe
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %in% c("Fall 2022")
)
mm_lefse_6 <- run_lefse(divexp_phy_small, group = "Cover", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

gps <- unique(mm_lefse_6@marker_table$enrich_group)

p_6 <- plot_cladogram(mm_lefse_6, color = covercols[gps], 
                      alpha = 0.5, clade_label_level = 5) + 
  theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none") 



### Spring 2023 LEfSe
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %in% c("Spring 2023")
)
mm_lefse_7 <- run_lefse(divexp_phy_small, group = "Cover", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

gps <- unique(mm_lefse_7@marker_table$enrich_group)

p_7 <- plot_cladogram(mm_lefse_7, color = covercols[gps], 
                      alpha = 0.5, clade_label_level = 5) + 
  theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none") 


### Summer 2023 LEfSe
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %in% c("Summer 2023")
)
mm_lefse_8 <- run_lefse(divexp_phy_small, group = "Cover", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

gps <- unique(mm_lefse_8@marker_table$enrich_group)

p_8 <- plot_cladogram(mm_lefse_8, color = covercols[gps], 
                      alpha = 0.5, clade_label_level = 5) + 
  theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none") 


### Fall 2023 LEfSe
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %in% c("Fall 2023")
)
mm_lefse_9 <- run_lefse(divexp_phy_small, group = "Cover", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

gps <- unique(mm_lefse_9@marker_table$enrich_group)

p_9 <- plot_cladogram(mm_lefse_9, color = covercols[gps], 
                      alpha = 0.5, clade_label_level = 5) + 
  theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none") 



### Spring 2024 LEfSe
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %in% c("Spring 2024")
)
mm_lefse_10 <- run_lefse(divexp_phy_small, group = "Cover", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

gps <- unique(mm_lefse_10@marker_table$enrich_group)

p_10 <- plot_cladogram(mm_lefse_10, color = covercols[gps], 
                      alpha = 0.5, clade_label_level = 5) + 
  theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none") 


### Summer 2024 LEfSe
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %in% c("Summer 2024")
)
mm_lefse_11 <- run_lefse(divexp_phy_small, group = "Cover", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

gps <- unique(mm_lefse_11@marker_table$enrich_group)

p_11 <- plot_cladogram(mm_lefse_11, color = covercols[gps], 
                      alpha = 0.5, clade_label_level = 5) + 
  theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none") 


### Fall 2024 LEfSe
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %in% c("Fall 2024")
)
mm_lefse_12 <- run_lefse(divexp_phy_small, group = "Cover", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

gps <- unique(mm_lefse_12@marker_table$enrich_group)

p_12 <- plot_cladogram(mm_lefse_12, color = covercols[gps], 
                      alpha = 0.5, clade_label_level = 5, ) + 
  theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none") 


### export plots for all timepoints

# p_all <- ggpubr::ggarrange(p_1, NA, NA,
#                   p_4, p_5, NA,
#                   NA, NA, NA,
#                   p_10, NA, NA, nrow = 4, ncol = 3, 
#                   labels = c("A) Spring 2021", "B) Summer 2021", "C) Fall 2021", 
#                              "D) Spring 2022", "E) Summer 2022", "F) Fall 2022", 
#                              "G) Spring 2023", "H) Summer 2023", "I) Fall 2023", 
#                              "J) Spring 2024", "K) Summer 2024", "L) Fall 2024"))
# ggpubr::ggexport(p_all, width = 6000, height = 6000, res = reso, 
#                  filename = paste0("Figures/microbial-diversity/lefse/", region, "_by season_cover.png"))


# export lefse data
write.csv(data.frame(mm_lefse_0@marker_table[, c(1:4)]), paste("Model-output/lefse/", region, "_overall_cover.csv"))
write.csv(data.frame(mm_lefse_1@marker_table[, c(1:4)]), paste("Model-output/lefse/", region, "_Spring 2021_cover.csv"))
write.csv(data.frame(mm_lefse_2@marker_table[, c(1:4)]), paste("Model-output/lefse/", region, "_Summer 2021_cover.csv"))
write.csv(data.frame(mm_lefse_3@marker_table[, c(1:4)]), paste("Model-output/lefse/", region, "_Fall 2021_cover.csv"))
write.csv(data.frame(mm_lefse_4@marker_table[, c(1:4)]), paste("Model-output/lefse/", region, "_Spring 2022_cover.csv"))
write.csv(data.frame(mm_lefse_5@marker_table[, c(1:4)]), paste("Model-output/lefse/", region, "_Summer 2022_cover.csv"))
write.csv(data.frame(mm_lefse_6@marker_table[, c(1:4)]), paste("Model-output/lefse/", region, "_Fall 2022_cover.csv"))
write.csv(data.frame(mm_lefse_7@marker_table[, c(1:4)]), paste("Model-output/lefse/", region, "_Spring 2023_cover.csv"))
write.csv(data.frame(mm_lefse_8@marker_table[, c(1:4)]), paste("Model-output/lefse/", region, "_Summer 2023_cover.csv"))
write.csv(data.frame(mm_lefse_9@marker_table[, c(1:4)]), paste("Model-output/lefse/", region, "_Fall 2023_cover.csv"))
write.csv(data.frame(mm_lefse_10@marker_table[, c(1:4)]), paste("Model-output/lefse/", region, "_Spring 2024_cover.csv"))
write.csv(data.frame(mm_lefse_11@marker_table[, c(1:4)]), paste("Model-output/lefse/", region, "_Summer 2024_cover.csv"))
write.csv(data.frame(mm_lefse_12@marker_table[, c(1:4)]), paste("Model-output/lefse/", region, "_Fall 2024_cover.csv"))
















### LEfSe with all experimental timepoints included
'%notin%' <- Negate('%in%')
# subset to timepoint
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %notin% c("Fall 2020")
)
# lefse analysis
mm_lefse_0 <- run_lefse(divexp_phy_small, group = "Cropping.system", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

# which treatments are represented?
gps <- unique(mm_lefse_0@marker_table$enrich_group)

# plot with cladogram
p_overall <- plot_cladogram(mm_lefse_0, color = cropsyscols[gps], 
                            alpha = 0.5, clade_label_level = 5) +
  theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none")  

# ggpubr::ggexport(p_overall + labs(tag = "A) Bacteria"), width = 4000, height = 3000, res = reso, 
#                  filename = paste0("Figures/microbial-diversity/lefse/", region, "_overall_Cropping.system.png"))





### Spring 2021 LEfSe
# subset to timepoint
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %in% c("Spring 2021")
)

# lefse analysis
mm_lefse_1 <- run_lefse(divexp_phy_small, group = "Cropping.system", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

# which treatments are represented?
gps <- unique(mm_lefse_1@marker_table$enrich_group)

# plot with cladogram
p_1 <- plot_cladogram(mm_lefse_1, color = cropsyscols[gps], 
                      alpha = 0.5, clade_label_level = 5) + 
  theme(plot.margin = margin(0, 0, 0, 0)) + 
  guides(fill = "none") #, shape = guide_legend(override.aes = list(ncol=1, color = cropsyscols[mm_lefse_1@marker_table$enrich_group][-1], shape = 15, size = 5)))    


### Summer 2021 LEfSe
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %in% c("Summer 2021")
)
mm_lefse_2 <- run_lefse(divexp_phy_small, group = "Cropping.system", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

gps <- unique(mm_lefse_2@marker_table$enrich_group)

p_2 <- plot_cladogram(mm_lefse_2, color = cropsyscols[gps], 
                      alpha = 0.5, clade_label_level = 5) + 
  theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none") 


### Fall 2021 LEfSe
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %in% c("Fall 2021")
)
mm_lefse_3 <- run_lefse(divexp_phy_small, group = "Cropping.system", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

gps <- unique(mm_lefse_3@marker_table$enrich_group)

p_3 <- plot_cladogram(mm_lefse_3, color = cropsyscols[gps], 
                      alpha = 0.5, clade_label_level = 5) + 
  theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none") 


### Spring 2022 LEfSe
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %in% c("Spring 2022")
)
mm_lefse_4 <- run_lefse(divexp_phy_small, group = "Cropping.system", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

gps <- unique(mm_lefse_4@marker_table$enrich_group)

p_4 <- plot_cladogram(mm_lefse_4, color = cropsyscols[gps], 
                      alpha = 0.5, clade_label_level = 5) + 
  theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none") 


### Summer 2022 LEfSe
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %in% c("Summer 2022")
)
mm_lefse_5 <- run_lefse(divexp_phy_small, group = "Cropping.system", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

gps <- unique(mm_lefse_5@marker_table$enrich_group)

p_5 <- plot_cladogram(mm_lefse_5, color = cropsyscols[gps], 
                      alpha = 0.5, clade_label_level = 5) + 
  theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none") 


### Fall 2022 LEfSe
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %in% c("Fall 2022")
)
mm_lefse_6 <- run_lefse(divexp_phy_small, group = "Cropping.system", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

gps <- unique(mm_lefse_6@marker_table$enrich_group)

p_6 <- plot_cladogram(mm_lefse_6, color = cropsyscols[gps], 
                      alpha = 0.5, clade_label_level = 5) + 
  theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none") 



### Spring 2023 LEfSe
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %in% c("Spring 2023")
)
mm_lefse_7 <- run_lefse(divexp_phy_small, group = "Cropping.system", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

gps <- unique(mm_lefse_7@marker_table$enrich_group)

p_7 <- plot_cladogram(mm_lefse_7, color = cropsyscols[gps], 
                      alpha = 0.5, clade_label_level = 5) + 
  theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none") 


### Summer 2023 LEfSe
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %in% c("Summer 2023")
)
mm_lefse_8 <- run_lefse(divexp_phy_small, group = "Cropping.system", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

gps <- unique(mm_lefse_8@marker_table$enrich_group)

p_8 <- plot_cladogram(mm_lefse_8, color = cropsyscols[gps], 
                      alpha = 0.5, clade_label_level = 5) + 
  theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none") 


### Fall 2023 LEfSe
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %in% c("Fall 2023")
)
mm_lefse_9 <- run_lefse(divexp_phy_small, group = "Cropping.system", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

gps <- unique(mm_lefse_9@marker_table$enrich_group)

p_9 <- plot_cladogram(mm_lefse_9, color = cropsyscols[gps], 
                      alpha = 0.5, clade_label_level = 5) + 
  theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none") 



### Spring 2024 LEfSe
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %in% c("Spring 2024")
)
mm_lefse_10 <- run_lefse(divexp_phy_small, group = "Cropping.system", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

gps <- unique(mm_lefse_10@marker_table$enrich_group)

p_10 <- plot_cladogram(mm_lefse_10, color = cropsyscols[gps], 
                       alpha = 0.5, clade_label_level = 5) + 
  theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none") 


### Summer 2024 LEfSe
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %in% c("Summer 2024")
)
mm_lefse_11 <- run_lefse(divexp_phy_small, group = "Cropping.system", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

gps <- unique(mm_lefse_11@marker_table$enrich_group)

p_11 <- plot_cladogram(mm_lefse_11, color = cropsyscols[gps], 
                       alpha = 0.5, clade_label_level = 5) + 
  theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none") 


### Fall 2024 LEfSe
divexp_phy_small <- phyloseq::subset_samples(
  divexp_phy,
  Time %in% c("Fall 2024")
)
mm_lefse_12 <- run_lefse(divexp_phy_small, group = "Cropping.system", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)

gps <- unique(mm_lefse_12@marker_table$enrich_group)

p_12 <- plot_cladogram(mm_lefse_12, color = cropsyscols[gps], 
                       alpha = 0.5, clade_label_level = 5, ) + 
  theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none") 


### export plots for all timepoints
p_all <- ggpubr::ggarrange(p_4, p_6, 
                           p_9,
                           nrow = 1, ncol = 3, 
                           labels = c("A) Spring 2022",  "B) Fall 2022", 
                                      "C) Fall 2023"))
ggpubr::ggexport(p_all, width = 6000, height = 2000, res = reso, 
                 filename = paste0("Figures/microbial-diversity/lefse/", region, "_by season_Cropping.system.png"))


# export lefse data
write.csv(data.frame(mm_lefse_0@marker_table[, c(1:4)]), paste("Model-output/lefse/", region, "_overall_Cropping.system.csv"))
write.csv(data.frame(mm_lefse_1@marker_table[, c(1:4)]), paste("Model-output/lefse/", region, "_Spring 2021_Cropping.system.csv"))
write.csv(data.frame(mm_lefse_2@marker_table[, c(1:4)]), paste("Model-output/lefse/", region, "_Summer 2021_Cropping.system.csv"))
write.csv(data.frame(mm_lefse_3@marker_table[, c(1:4)]), paste("Model-output/lefse/", region, "_Fall 2021_Cropping.system.csv"))
write.csv(data.frame(mm_lefse_4@marker_table[, c(1:4)]), paste("Model-output/lefse/", region, "_Spring 2022_Cropping.system.csv"))
write.csv(data.frame(mm_lefse_5@marker_table[, c(1:4)]), paste("Model-output/lefse/", region, "_Summer 2022_Cropping.system.csv"))
write.csv(data.frame(mm_lefse_6@marker_table[, c(1:4)]), paste("Model-output/lefse/", region, "_Fall 2022_Cropping.system.csv"))
write.csv(data.frame(mm_lefse_7@marker_table[, c(1:4)]), paste("Model-output/lefse/", region, "_Spring 2023_Cropping.system.csv"))
write.csv(data.frame(mm_lefse_8@marker_table[, c(1:4)]), paste("Model-output/lefse/", region, "_Summer 2023_Cropping.system.csv"))
write.csv(data.frame(mm_lefse_9@marker_table[, c(1:4)]), paste("Model-output/lefse/", region, "_Fall 2023_Cropping.system.csv"))
write.csv(data.frame(mm_lefse_10@marker_table[, c(1:4)]), paste("Model-output/lefse/", region, "_Spring 2024_Cropping.system.csv"))
write.csv(data.frame(mm_lefse_11@marker_table[, c(1:4)]), paste("Model-output/lefse/", region, "_Summer 2024_Cropping.system.csv"))
write.csv(data.frame(mm_lefse_12@marker_table[, c(1:4)]), paste("Model-output/lefse/", region, "_Fall 2024_Cropping.system.csv"))








