
library(tidyverse)
library(pheatmap)




#### diversity anova
rm(list=ls())
regions <- c("16S", "ITS", "AMF", "plant pathogen", "saprotroph")

# import anova data for each microbial group
resp_1 <- read.csv(paste0("Model-output/anova_", regions[1], "/all-anova.csv"))
resp_1$name <- rep(regions[1], dim(resp_1)[1])
resp_2 <- read.csv(paste0("Model-output/anova_", regions[2], "/all-anova.csv"))
resp_2$name <- rep(regions[2], dim(resp_2)[1])
resp_3 <- read.csv(paste0("Model-output/anova_", regions[3], "/all-anova.csv"))
resp_3$name <- rep(regions[3], dim(resp_3)[1])
resp_4 <- read.csv(paste0("Model-output/anova_", regions[4], "/all-anova.csv"))
resp_4$name <- rep(regions[4], dim(resp_4)[1])
resp_5 <- read.csv(paste0("Model-output/anova_", regions[5], "/all-anova.csv"))
resp_5$name <- rep(regions[5], dim(resp_5)[1])

# combine all data
resp <- rbind(resp_1, resp_2, resp_3,  resp_4,  resp_5)

# filter for richness and select anova results (asterisks for significance)
resp_f <- resp %>% 
  filter(resp %in% c("bac_richness", "fun_richness", 
                     "amf_richness", "amf_relabund", 
                     "plantpath_richness","plantpath_relabund",
                     "sap_richness", "sap_relabund")) %>% 
  dplyr::select(timepoint, resp, cover_sig, cropsys_sig, cover.cropsys_sig, cover, cropsys, cover.cropsys)

# paste together f-value and significance
resp_f <- resp_f %>% 
  mutate(cover_f = stringr::str_split(cover, ", ") %>% map_chr(., 1)) %>% 
  mutate(cover_f = stringr::str_split(cover_f, "F=") %>% map_chr(., 2)) %>% 
  mutate(cropsys_f = stringr::str_split(cropsys, ", ") %>% map_chr(., 1)) %>% 
  mutate(cropsys_f = stringr::str_split(cropsys_f, "F=") %>% map_chr(., 2)) %>% 
  mutate(cover.cropsys_f = stringr::str_split(cover.cropsys, ", ") %>% map_chr(., 1)) %>% 
  mutate(cover.cropsys_f = stringr::str_split(cover.cropsys_f, "F=") %>% map_chr(., 2)) %>% 
  mutate(cover_effect = paste(cover_f, cover_sig) ) %>% 
  mutate(cropsys_effect = paste(cropsys_f, cropsys_sig) ) %>% 
  mutate(cover.cropsys_effect = paste(cover.cropsys_f, cover.cropsys_sig) )

# # remove fall 2020
# resp_f <- resp_f[-which(resp_f$timepoint=="Fall 2020"),]


# cbind different microbial groups
resp_h <- (resp_f) %>% 
  dplyr::select(resp, timepoint, cover_effect, cropsys_effect, cover.cropsys_effect) %>% 
  pivot_wider(names_from=resp, values_from=c(cover_effect, cropsys_effect, cover.cropsys_effect))

# reorder columns
resp_e <- resp_h[,c(1,
            which(str_detect(colnames(resp_h), pattern="bac")), 
            which(str_detect(colnames(resp_h), pattern="fun")), 
            which(str_detect(colnames(resp_h), pattern="amf_richness")), 
            which(str_detect(colnames(resp_h), pattern="amf_relabund")), 
            which(str_detect(colnames(resp_h), pattern="plantpath_richness")), 
            which(str_detect(colnames(resp_h), pattern="plantpath_relabund")), 
            which(str_detect(colnames(resp_h), pattern="sap_richness")),
            which(str_detect(colnames(resp_h), pattern="sap_relabund")))]

resp_t <- t(resp_e)


write.csv(resp_t, "Model-output/*tables/diversity_sig.csv")




#### soil anova
rm(list=ls())
props <- c("GMC", "noN", "nhN",  "WEC", "WEN", "MBC", "BG", "NAG", "PHOS") # removed "inorgN",

resp_1 <- read.csv(paste0("Model-output/anova_16S/all-anova.csv"))

resp <- rbind(resp_1)

# filter for richness and select anova results
resp_f <- resp %>% 
  filter(resp %in% props) %>% 
  dplyr::select(timepoint, resp, cover, cropsys, cover.cropsys, rsq)

# remove fall 2020
resp_f <- resp_f[-1,]


write.csv(t(resp_f), "Model-output/*tables/soil.csv")

# filter for richness and select anova results (asterisks for significance)
resp_f <- resp %>% 
  filter(resp %in% props) %>% 
  dplyr::select(timepoint, resp, cover_sig, cropsys_sig, cover.cropsys_sig, cover, cropsys, cover.cropsys)

# paste together f-value and significance
resp_s <- resp_f %>% 
  mutate(cover_f = stringr::str_split(cover, ", ") %>% map_chr(., 1)) %>% 
  mutate(cover_f = stringr::str_split(cover_f, "F=") %>% map_chr(., 2)) %>% 
  mutate(cropsys_f = stringr::str_split(cropsys, ", ") %>% map_chr(., 1)) %>% 
  mutate(cropsys_f = stringr::str_split(cropsys_f, "F=") %>% map_chr(., 2)) %>% 
  mutate(cover.cropsys_f = stringr::str_split(cover.cropsys, ", ") %>% map_chr(., 1)) %>% 
  mutate(cover.cropsys_f = stringr::str_split(cover.cropsys_f, "F=") %>% map_chr(., 2)) %>% 
  mutate(cover_effect = paste(cover_f, cover_sig) ) %>% 
  mutate(cropsys_effect = paste(cropsys_f, cropsys_sig) ) %>% 
  mutate(cover.cropsys_effect = paste(cover.cropsys_f, cover.cropsys_sig) )

# remove fall 2020
resp_f <- resp_s[-which(str_detect(resp_s$timepoint, pattern = "Fall 2020")),]

# change to short format
resp_f_s <- resp_f  %>% 
  dplyr::select(resp, timepoint, cover_effect, cropsys_effect, cover.cropsys_effect) %>% 
  pivot_wider(names_from = resp, 
              values_from = c(cover_effect, cropsys_effect, cover.cropsys_effect))

# transpose
resp_f_s <- t(resp_f_s) %>% 
          janitor::row_to_names(row_number=1)

# reorder
resp_f_s <- as.data.frame(resp_f_s)
resp_f_s$nm <- rownames(resp_f_s)
resp_f_s <- resp_f_s %>% 
            mutate(response = stringr::str_split(nm, "_") %>% 
                     map_chr(., 3))
resp_f_s <- resp_f_s %>% 
  arrange(response) 


write.csv(resp_f_s, "Model-output/*tables/soil_sig.csv")





#### permanova 
rm(list=ls())
regions <- c("Bacteria", "Fungi", "Fungi_AMF", "Fungi_plant pathogen", "Fungi_saprotroph")

# read in data for each microbial group
resp_1 <- read.csv(paste0("Model-output/permanova/", regions[1], "_permanova_by season_*.csv"), row.names=2)
resp_2 <- read.csv(paste0("Model-output/permanova/", regions[2], "_permanova_by season_*.csv"), row.names=2)
resp_3 <- read.csv(paste0("Model-output/permanova/", regions[3], "_permanova_by season_*.csv"), row.names=2)
resp_4 <- read.csv(paste0("Model-output/permanova/", regions[4], "_permanova_by season_*.csv"), row.names=2)
resp_5 <- read.csv(paste0("Model-output/permanova/", regions[5], "_permanova_by season_*.csv"), row.names=2)

# combine data
resp <- rbind(resp_1, resp_2, resp_3, resp_4, resp_5)


# remove fall 2020
resp_f <- resp[,-c(1:2)]

# transpose
resp_t <- as.data.frame(t(resp_f))

# select asterisks only 
resp_f_s <- resp_t %>% 
  mutate(across(everything(),~ gsub(" ","% ", .))) %>% 
  mutate(across(everything(),~ gsub("NA","", .))) 

# # select last three columns
# resp_f_s <- as.data.frame(resp_f_s)
# resps <- resp_f_s[, c((dim(resp_f_s)[2]/2+1):dim(resp_f_s)[2])]

write.csv(t(resp_f_s), paste0("Model-output/*tables/permanova_sig.csv"))





#### dbrda: variance by predictors
rm(list=ls())
regions <- c("Bacteria", "Fungi", "Fungi_AMF", "Fungi_plant pathogen", "Fungi_saprotroph")
labels_main <- c("A) Bacteria", " B) Fungi", "C) Arbuscular mycorrhizal fungi", "D) Plant pathogenic fungi", "E) Saprotrophic fungi")
legends <- c(F, T, F, F, F)
widths <- c(1600, 1700, 1600, 1600, 1600)
seasons <- c("Fall 2020", 
             "Spring 2021", "Summer 2021", "Fall 2021",
             "Spring 2022", "Summer 2022", "Fall 2022",
             "Spring 2023", "Summer 2023", "Fall 2023",
             "Spring 2024", "Summer 2024", "Fall 2024")



for(m in 1:length(regions)){
  # read in data for each season
  resps <- list(NA)
  for(i in 1:length(seasons)){
    resp_1 <- read.csv(paste0("Model-output/db-rda/", regions[m], "_all_variance-by-predictors_", seasons[i], ".csv"))
    resps[[i]] <- resp_1
  }
  resps
  resps3 <- bind_rows(resps, .id = "timepoint")
  
  # extract r2 for each soil property
  resps4 <- resps3[,c("timepoint", "X", "r2")]
  resps5 <- spread(data=resps4, key=timepoint, value=c(r2))
  resps6 <- resps5[,-1]
  rownames(resps6) <- resps5$X
  resps6 <- resps6[, c(1, 6:13, 2:5)]
  resp_r <- as.data.frame(t(resps6))
  rownames(resp_r) <- seasons
  
  # extract pvalue for each soil property 
  resps4 <- resps3[,c("timepoint", "X", "pvalue")]
  resps5 <- spread(data=resps4, key=timepoint, value=c(pvalue))
  resps6 <- resps5[,-1]
  rownames(resps6) <- resps5$X
  resps6 <- resps6[, c(1, 6:13, 2:5)]
  resp_p <- as.data.frame(t(resps6))
  rownames(resp_p) <- seasons
  # replace pvalues with asterisks
  resp_s <- resp_p
  makeStars <- function(x){
    stars <- c("****", "***", "**", "*", "")
    vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1.01)
    i <- findInterval(x, vec)
    stars[i]
  }
  for(k in 1:dim(resp_s)[2]){
    resp_s[,k] <- makeStars(resp_s[,k])
  }

  
  p <- pheatmap(as.matrix(t(resp_r)), display_numbers = as.matrix(t(resp_s)), legend = legends[m], 
                color = colorRampPalette(c("white", "cyan4"))(100), legend.breaks = seq(0, 1, by = 0.25), legend.labels = seq(0, 1, by = 0.25),
                cluster_rows = F, cluster_cols = F, main=labels_main[m])
  
  save_pheatmap_png <- function(x, filename, width=widths, height=1000, res = 150) {
    png(filename, width = width, height = height, res = res)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  
  save_pheatmap_png(p, filename = paste0("Figures/microbial-diversity/rda/", regions[m], ".png"), width=widths[m], height=1200, res=300)
  
  write.csv(t(resp_r), paste0("Model-output/*tables/dbrda_variance-by-predictors_", regions[m], ".csv"))
}  







#### dbrda: total variance
rm(list=ls())
regions <- c("Bacteria", "Fungi", "Fungi_AMF", "Fungi_plant pathogen", "Fungi_saprotroph")
seasons <- c("Fall 2020", 
             "Spring 2021", "Summer 2021", "Fall 2021",
             "Spring 2022", "Summer 2022", "Fall 2022",
             "Spring 2023", "Summer 2023", "Fall 2023",
             "Spring 2024", "Summer 2024", "Fall 2024")

sigdf <- data.frame(season = seasons,
                    Bacteria = rep(NA, length(seasons)),
                    Fungi = rep(NA, length(seasons)),
                    Fungi_AMF = rep(NA, length(seasons)),
                    'Fungi_plant pathogen' = rep(NA, length(seasons)),
                    Fungi_saprotroph = rep(NA, length(seasons)))

for(i in 1:length(regions)){
  sigs <- read.csv(paste0("Model-output/db-rda/*",regions[i],"_dbrda_by timepoint_total-r2-sig.csv"), row.names = 1)
  df_rda <- read.csv(paste0("Model-output/db-rda/*", regions[i], "_dbrda_by timepoint.csv"), row.names = 1)
  colnames(df_rda) <- df_rda[1,]
  df_rda <- df_rda[-1,]
  df_rda <- df_rda %>% 
    mutate(model.r = stringr::str_split(model, pattern="R2=") %>%  map_chr(., 2)) 
  df_rda$model.r <- as.numeric(df_rda$model.r)
  df_rda$sig <- c(sigs[1,c(1:length(seasons))])
  sigdf[, i+1] <- paste0(df_rda$model.r*100, "% ", df_rda$sig)
}


write.csv(sigdf, paste0("Model-output/*tables/dbrda_by timepoint_total-r2-sig.csv"))
