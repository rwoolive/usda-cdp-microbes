
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("lefser")

browseVignettes("lefser")


load("Raw-data/sequence/All_agmicrobiome/Fungilist.RData")
