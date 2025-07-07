







library(vegan)
library(tidyverse)
library(ecodist)
library(ape)
library(lme4)
library(emmeans)
library(nlme)
library(multcomp)





##################### Import soil data #####################
# field plot soil property data
dat <- read.csv("Processed-data/WTREC-CDP_processed-data-simple.csv")
# subset dat to 0-10cm depth
dat <- dat[which(dat$Depth=="0-10 cm"),] 
dat$Season.Year <- paste(dat$Season, dat$Year, sep=" ")
dat$Season.Year <- as.factor(dat$Season.Year )
dat$Season.Year <- factor(dat$Season.Year, levels(dat$Season.Year)[c(1,6,10,2,7,11,3,8,12,4,9,13,5)])


# make factors
dat <- dat %>% mutate_at(c("Rep", "Year", "Season", "Cropping.system", "Cover"), factor)
dat$Cropping.system <- factor(dat$Cropping.system, levels(dat$Cropping.system)[c(1,4,3,2)])
dat$Cover <- factor(dat$Cover, levels(dat$Cover)[c(2,4,1,5,3)])

################## sampling dates
nsamps <- 13
dat %>%
  dplyr::group_by(Season.Year) %>%
  dplyr::summarize(mean=mean(GMC, na.rm=T))
################## 

cropsyscols <- c("darkgoldenrod4", "darkolivegreen", "darkolivegreen3", "deepskyblue3")
covercols <- c("chocolate", "goldenrod", "forestgreen", "darkcyan", "darkorchid4")
ssnyrcols <- c("white", 
               RColorBrewer::brewer.pal(n = 4, name = "Purples")[2:4], 
               RColorBrewer::brewer.pal(n = 4, name = "YlGn")[2:4], 
               RColorBrewer::brewer.pal(n = 4, name = "Reds")[2:4], 
               RColorBrewer::brewer.pal(n = 4, name = "Blues")[2:4])
yrcols <- c("white", "purple", "aquamarine3", "brown2", "blue")





# response variables:
dat <- dat %>% mutate(across(c("GMC", "noN", "nhN", "inorgN", "WEC", "WEN", "MBC", "BG", "NAG", "PHOS"), as.numeric))

# select data
dat <- dat %>% 
  dplyr::select(Rep, Depth, Year, Season.Year, Season, Cropping.system, Cover, GMC, noN, nhN, inorgN, WEC, WEN, MBC, BG, NAG, PHOS)

dat_cover <- dat %>% 
  group_by(Cover, Season.Year) %>% 
  dplyr::summarise_each(funs(mean(., na.rm=T), se=sd(., na.rm=T)/sqrt(n())), GMC:PHOS)

dat_cropsys <- dat %>% 
  group_by(Cropping.system, Season.Year) %>% 
  dplyr::summarise_each(funs(mean(., na.rm=T), se=sd(., na.rm=T)/sqrt(n())), GMC:PHOS)

dat_cover_cropsys <- dat %>% 
  group_by(Cover, Cropping.system, Season.Year) %>% 
  dplyr::summarise_each(funs(mean(., na.rm=T), se=sd(., na.rm=T)/sqrt(n())), GMC:PHOS)




seasons <- c("Fall 2020", 
             "Spring 2021", "Summer 2021", "Fall 2021",
             "Spring 2022", "Summer 2022", "Fall 2022",
             "Spring 2023", "Summer 2023", "Fall 2023",
             "Spring 2024", "Summer 2024", "Fall 2024")

resps <- c("GMC", "noN", "nhN", "inorgN", "WEC", "WEN", "MBC", "BG", "NAG", "PHOS")
lets <- LETTERS[(1:length(resps))]
mains <- c(expression(paste("Gravimetric moisture content (g g"^"-1",")")), 
           expression(paste("Nitrate-N (mg kg"^-1,")")),
           expression(paste("Ammonium-N (mg kg"^-1,")")),
           expression(paste("Inorganic N (mg kg"^-1,")")),
           expression(paste("Water-extractable OC (mg kg"^-1,")")),
           expression(paste("Water-extractable N (mg kg"^-1,")")),
           expression(paste("Microbial biomass C (mg kg"^-1,")")), 
           expression(paste(beta,"-glucosidase activity (nmol g"^-1," h"^-1,")")),
           expression(paste("N-acetyl-", beta, "-glucosaminidase activity (nmol g"^-1," h"^-1,")")),
           expression(paste("Phosphatase activity (nmol g"^-1," h"^-1,")")))
regions <- rep("16S", length(resps))
heights <- c(rep(2.5, 9), 2.5)
nox <- rep("", length(seasons))
xaxes <- list(nox, nox, nox, nox,
              nox, nox, nox, nox, nox, nox)
widths <- c(8.5, rep(8.5,9))
legs <- c( T, rep(F, 9))



#######################################
# figures of microbial diversity across seasons, colored by cover
#######################################

for(i in 1:length(resps)){
  ## import anova data ## 
df <- read.csv(paste0("Model-output/anova_", regions[i], "/all-anova.csv"), row.names=1)
# filter to specific response
df <- df %>% filter(resp == resps[i])
# determine which timepoints are significant
sigs <- which(df$cover_sig %in% c("*", "**", "***"))
# select data for plot
dat_sub <- dat_cover %>% dplyr::select("Season.Year", "Cover", paste0(resps[i], "_mean"), paste0(resps[i], "_se")) 
colnames(dat_sub) <- c("Season.Year", "Cover", "Response", "SE")
# calculate where to put the asterisks
dat_max <- dat_sub %>% 
  group_by(Season.Year) %>% 
  dplyr::summarize(max = max(Response, na.rm=T))

# plot
p <- ggplot(dat_sub, aes(x=as.numeric(Season.Year), y=Response, color=Cover)) +
  geom_vline(xintercept = c(2,5,8,11), color="gray77", linetype="dotted") +
  #annotate(geom = "text", x = sigs, y = dat_max$max[sigs]+0.25*(max(dat_max$max)-min(dat_max$max)), label = df$cover_sig[sigs], size=5) +
  geom_point() +
  geom_line() + 
  scale_color_manual(values=covercols) +
  labs(color="Cover", x="", y="", title=mains[i]) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title.y = element_text(size=8)) +
  scale_x_continuous(breaks=c(1:13), labels = xaxes[[i]]) 
# if(legs[i]==F) {q <- p + guides(color=F)
#             # export plot
#             ggpubr::ggexport(q, height=heights[i], width=widths[i], 
#                              filename = paste0("Figures/microbial-diversity/diversity/Cover ", lets[i], ") ", resps[i], ".pdf"))
# } 
# else(# export plot
  ggpubr::ggexport(p, height=heights[i], width=widths[i], 
                   filename = paste0("Figures/microbial-diversity/diversity/Cover ", lets[i], ") ", resps[i], ".pdf"))#)
}





#######################################
# figures of microbial diversity across seasons, colored by cropping.system
#######################################

for(i in 1:length(resps)){
  ## import anova data ## 
  df <- read.csv(paste0("Model-output/anova_", regions[i], "/all-anova.csv"), row.names=1)
  # filter to specific response
  df <- df %>% filter(resp == resps[i])
  # determine which timepoints are significant
  sigs <- which(df$cropsys_sig %in% c("*", "**", "***"))
  # select data for plot
  dat_sub <- dat_cropsys %>% dplyr::select("Season.Year", "Cropping.system", paste0(resps[i], "_mean"), paste0(resps[i], "_se")) 
  colnames(dat_sub) <- c("Season.Year", "Cropping.system", "Response", "SE")
  # calculate where to put the asterisks
  dat_max <- dat_sub %>% 
    group_by(Season.Year) %>% 
    dplyr::summarize(max = max(Response, na.rm=T))
  
  # plot
  p <- ggplot(dat_sub, aes(x=as.numeric(Season.Year), y=Response, color=Cropping.system)) +
    geom_vline(xintercept = c(2,5,8,11), color="gray77", linetype="dotted") +
    #annotate(geom = "text", x = sigs, y = dat_max$max[sigs]+0.25*(max(dat_max$max)-min(dat_max$max)), label = df$cropsys_sig[sigs], size=5) +
    geom_point() +
    geom_line() + 
    scale_color_manual(values=cropsyscols) +
    labs(color="Cropping system", x="", y="", title=mains[i]) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title.y = element_text(size=8)) +
    scale_x_continuous(breaks=c(1:13), labels = xaxes[[i]]) 
  # if(legs[i]==F) {q <- p + guides(color=F)
  # # export plot
  # ggpubr::ggexport(q, height=heights[i], width=widths[i], 
  #                  filename = paste0("Figures/microbial-diversity/diversity/Cropping system ", lets[i], ") ", resps[i], ".pdf"))
  # } 
  # else(# export plot
    ggpubr::ggexport(p, height=heights[i], width=widths[i], 
                     filename = paste0("Figures/microbial-diversity/diversity/Cropping system ", lets[i], ") ", resps[i], ".pdf"))#)
}





#######################################
# figures of microbial diversity across seasons, colored by cover and faceted by cropsys
#######################################

for(i in 1:length(resps)){
  ## import anova data ## 
  df <- read.csv(paste0("Model-output/anova_", regions[i], "/all-anova.csv"), row.names=1)
  # filter to specific response
  df <- df %>% filter(resp == resps[i])
  # determine which timepoints are significant
  sigs <- which(df$cover.cropsys_sig %in% c("*", "**", "***"))
  # select data for plot
  dat_sub <- dat_cover_cropsys %>% dplyr::select("Season.Year", "Cover", "Cropping.system", paste0(resps[i], "_mean"), paste0(resps[i], "_se")) 
  colnames(dat_sub) <- c("Season.Year", "Cover", "Cropping.system", "Response", "SE")
  # calculate where to put the asterisks
  dat_max <- dat_sub %>% 
    group_by(Season.Year) %>% 
    dplyr::summarize(max = max(Response, na.rm=T))
  
  # plot
  p <- ggplot(dat_sub, aes(x=as.numeric(Season.Year), y=Response, color=Cover)) +
    geom_vline(xintercept = c(2,5,8,11), color="gray77", linetype="dotted") +
    #annotate(geom = "text", x = sigs, y = dat_max$max[sigs]+0.25*(max(dat_max$max)-min(dat_max$max)), label = df$cover_sig[sigs], size=5) +
    geom_point() +
    geom_line() + 
    scale_color_manual(values=covercols) +
    labs(color="Cover", x="", y="", title=mains[i]) +
    theme_bw() +
    facet_wrap(~Cropping.system, nrow=4) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title.y = element_text(size=8)) +
    scale_x_continuous(breaks=c(1:13), labels = xaxes[[8]]) 
# export plot
    ggpubr::ggexport(p, height=heights[8]*2.5, width=widths[1]*0.8, 
                     filename = paste0("Figures/microbial-diversity/diversity/CoverXCropsys ", lets[i], ") ", resps[i], ".pdf"))
}







# make factors
dat <- dat %>% mutate_at(c("Replicate", "Year", "Season", "Cropping.system", "Cover", "Season.Year"), factor)
dat$Cropping.system <- factor(dat$Cropping.system, levels(dat$Cropping.system)[c(1,4,3,2)])
dat$Cover <- factor(dat$Cover, levels(dat$Cover)[c(2,4,1,5,3)])
dat$Season.Year <- factor(dat$Season.Year, levels(dat$Season.Year)[c(1,6,10,2,7,11,3,8,12,4,9,13,5)])

################## sampling dates
nsamps <- 13
dat %>%
  dplyr::group_by(Season.Year) %>%
  dplyr::summarize(mean=mean(fun_richness, na.rm=T))
################## 

cropsyscols <- c("darkgoldenrod4", "darkolivegreen", "darkolivegreen3", "deepskyblue3")
covercols <- c("chocolate", "goldenrod", "forestgreen", "darkcyan", "darkorchid4")
ssnyrcols <- c("white", 
               RColorBrewer::brewer.pal(n = 4, name = "Purples")[2:4], 
               RColorBrewer::brewer.pal(n = 4, name = "YlGn")[2:4], 
               RColorBrewer::brewer.pal(n = 4, name = "Reds")[2:4], 
               RColorBrewer::brewer.pal(n = 4, name = "Blues")[2:4])
yrcols <- c("white", "purple", "aquamarine3", "brown2", "blue")





# response variables:
dat <- dat %>% mutate(across(c("bac_richness", #"bac_shannon.div", "bac_simpson.div", "bac_invsimpson.div", "bac_evenness",
                               "fun_richness", #"fun_shannon.div", "fun_simpson.div", "fun_invsimpson.div", "fun_evenness",
                               "amf_relabund", "amf_richness", #"amf_shannon.div", "amf_simpson.div", "amf_invsimpson.div", "amf_evenness",
                               "plantpath_relabund", "plantpath_richness", #"plantpath_shannon.div", "plantpath_simpson.div", "plantpath_invsimpson.div", "plantpath_evenness",
                               "sap_relabund", "sap_richness"), as.numeric))#, "sap_shannon.div", "sap_simpson.div", "sap_invsimpson.div", "sap_evenness"))


dat_cover <- dat %>% 
  group_by(Cover, Season.Year) %>% 
  dplyr::summarise_each(funs(mean(., na.rm=T), se=sd(., na.rm=T)/sqrt(n())), bac_richness:sap_relabund)

dat_cropsys <- dat %>% 
  group_by(Cropping.system, Season.Year) %>% 
  dplyr::summarise_each(funs(mean(., na.rm=T), se=sd(., na.rm=T)/sqrt(n())), bac_richness:sap_relabund)

dat_cover_cropsys <- dat %>% 
  group_by(Cover, Cropping.system, Season.Year) %>% 
  dplyr::summarise_each(funs(mean(., na.rm=T), se=sd(., na.rm=T)/sqrt(n())), bac_richness:sap_relabund)




seasons <- c("Fall 2020", 
             "Spring 2021", "Summer 2021", "Fall 2021",
             "Spring 2022", "Summer 2022", "Fall 2022",
             "Spring 2023", "Summer 2023", "Fall 2023",
             "Spring 2024", "Summer 2024", "Fall 2024")

resps <- c("bac_richness", #"bac_shannon.div", "bac_simpson.div", "bac_invsimpson.div", "bac_evenness",
           "fun_richness", #"fun_shannon.div", "fun_simpson.div", "fun_invsimpson.div", "fun_evenness",
           "amf_richness","amf_relabund",  #"amf_shannon.div", "amf_simpson.div", "amf_invsimpson.div", "amf_evenness",
           "plantpath_richness","plantpath_relabund",  #"plantpath_shannon.div", "plantpath_simpson.div", "plantpath_invsimpson.div", "plantpath_evenness",
           "sap_richness", "sap_relabund")
lets <- LETTERS[(1:length(resps))]
mains <- c("Bacteria richness", "Fungi richness", 
           "AMF richness", "AMF relative abundance", 
           "Plant pathogen richness", "Plant pathogen relative abundance", 
           "Saprotroph richness", "Saprotroph relative abundance" )
regions <- c("16S", "ITS", "AMF", "ITS", "plant pathogen", "ITS", "saprotroph", "ITS")
heights <- c(rep(2.5, 7), 3)
nox <- rep("", length(seasons))
xaxes <- list(nox, nox, nox, nox,
              nox, nox, nox, seasons)
widths <- c(8.5, rep(7,7))
legs <- c( T, rep(F, 7))



#######################################
# figures of microbial diversity across seasons, colored by cover
#######################################

for(i in 1:length(resps)){
  ## import anova data ## 
  df <- read.csv(paste0("Model-output/anova_", regions[i], "/all-anova.csv"), row.names=1)
  # filter to specific response
  df <- df %>% filter(resp == resps[i])
  # determine which timepoints are significant
  sigs <- which(df$cover_sig %in% c("*", "**", "***"))
  # select data for plot
  dat_sub <- dat_cover %>% dplyr::select("Season.Year", "Cover", paste0(resps[i], "_mean"), paste0(resps[i], "_se")) 
  colnames(dat_sub) <- c("Season.Year", "Cover", "Response", "SE")
  # calculate where to put the asterisks
  dat_max <- dat_sub %>% 
    group_by(Season.Year) %>% 
    dplyr::summarize(max = max(Response, na.rm=T))
  
  # plot
  p <- ggplot(dat_sub, aes(x=as.numeric(Season.Year), y=Response, color=Cover)) +
    geom_vline(xintercept = c(2,5,8,11), color="gray77", linetype="dotted") +
    annotate(geom = "text", x = sigs, y = dat_max$max[sigs]+0.25*(max(dat_max$max)-min(dat_max$max)), label = df$cover_sig[sigs], size=5) +
    geom_point() +
    geom_line() + 
    scale_color_manual(values=covercols) +
    labs(color="Cover", x="", y="", title=paste0(lets[i], ") ", mains[i])) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title.y = element_text(size=8)) +
    scale_x_continuous(breaks=c(1:13), labels = xaxes[[i]]) 
  if(legs[i]==F) {q <- p + guides(color=F)
  # export plot
  ggpubr::ggexport(q, height=heights[i], width=widths[i], 
                   filename = paste0("Figures/microbial-diversity/diversity/Cover ", lets[i], ") ", resps[i], ".pdf"))
  } 
  else(# export plot
    ggpubr::ggexport(p, height=heights[i], width=widths[i], 
                     filename = paste0("Figures/microbial-diversity/diversity/Cover ", lets[i], ") ", resps[i], ".pdf")))
}





#######################################
# figures of microbial diversity across seasons, colored by cropping.system
#######################################

for(i in 1:length(resps)){
  ## import anova data ## 
  df <- read.csv(paste0("Model-output/anova_", regions[i], "/all-anova.csv"), row.names=1)
  # filter to specific response
  df <- df %>% filter(resp == resps[i])
  # determine which timepoints are significant
  sigs <- which(df$cropsys_sig %in% c("*", "**", "***"))
  # select data for plot
  dat_sub <- dat_cropsys %>% dplyr::select("Season.Year", "Cropping.system", paste0(resps[i], "_mean"), paste0(resps[i], "_se")) 
  colnames(dat_sub) <- c("Season.Year", "Cropping.system", "Response", "SE")
  # calculate where to put the asterisks
  dat_max <- dat_sub %>% 
    group_by(Season.Year) %>% 
    dplyr::summarize(max = max(Response, na.rm=T))
  
  # plot
  p <- ggplot(dat_sub, aes(x=as.numeric(Season.Year), y=Response, color=Cropping.system)) +
    geom_vline(xintercept = c(2,5,8,11), color="gray77", linetype="dotted") +
    annotate(geom = "text", x = sigs, y = dat_max$max[sigs]+0.25*(max(dat_max$max)-min(dat_max$max)), label = df$cropsys_sig[sigs], size=5) +
    geom_point() +
    geom_line() + 
    scale_color_manual(values=cropsyscols) +
    labs(color="Cropping system", x="", y="", title=paste0(lets[i], ") ", mains[i])) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title.y = element_text(size=8)) +
    scale_x_continuous(breaks=c(1:13), labels = xaxes[[i]]) 
  if(legs[i]==F) {q <- p + guides(color=F)
  # export plot
  ggpubr::ggexport(q, height=heights[i], width=widths[i], 
                   filename = paste0("Figures/microbial-diversity/diversity/Cropping system ", lets[i], ") ", resps[i], ".pdf"))
  } 
  else(# export plot
    ggpubr::ggexport(p, height=heights[i], width=widths[i], 
                     filename = paste0("Figures/microbial-diversity/diversity/Cropping system ", lets[i], ") ", resps[i], ".pdf")))
}






#######################################
# figures of microbial diversity across seasons, colored by cover and faceted by cropsys
#######################################

for(i in 1:length(resps)){
  ## import anova data ## 
  df <- read.csv(paste0("Model-output/anova_", regions[i], "/all-anova.csv"), row.names=1)
  # filter to specific response
  df <- df %>% filter(resp == resps[i])
  # determine which timepoints are significant
  sigs <- which(df$cover.cropsys_sig %in% c("*", "**", "***"))
  # select data for plot
  dat_sub <- dat_cover_cropsys %>% dplyr::select("Season.Year", "Cover", "Cropping.system", paste0(resps[i], "_mean"), paste0(resps[i], "_se")) 
  colnames(dat_sub) <- c("Season.Year", "Cover", "Cropping.system", "Response", "SE")
  # calculate where to put the asterisks
  dat_max <- dat_sub %>% 
    group_by(Season.Year) %>% 
    dplyr::summarize(max = max(Response, na.rm=T))
  
  # plot
  p <- ggplot(dat_sub, aes(x=as.numeric(Season.Year), y=Response, color=Cover)) +
    geom_vline(xintercept = c(2,5,8,11), color="gray77", linetype="dotted") +
    #annotate(geom = "text", x = sigs, y = dat_max$max[sigs]+0.25*(max(dat_max$max)-min(dat_max$max)), label = df$cover_sig[sigs], size=5) +
    geom_point() +
    geom_line() + 
    scale_color_manual(values=covercols) +
    labs(color="Cover", x="", y="", title=paste0(mains[i])) +
    theme_bw() +
    facet_wrap(~Cropping.system, nrow=4) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title.y = element_text(size=8)) +
    scale_x_continuous(breaks=c(1:13), labels = xaxes[[8]]) 
  # export plot
  ggpubr::ggexport(p, height=heights[8]*2.5, width=widths[1]*0.8, 
                   filename = paste0("Figures/microbial-diversity/diversity/CoverXCropsys ", lets[i], ") ", resps[i], ".pdf"))
}






