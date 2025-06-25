##################################################
###### 2. ordination




library(vegan)
library(tidyverse)
library(ecodist)
library(ape)
library(lme4)
library(emmeans)
library(nlme)
library(multcomp)
library(ggpubr)
library(dplyr)
library(stringr)



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

# # insert POC, sand, silt, and clay values into soildat
# texture <- read.csv("Processed-data/WTREC-CDP_processedpooled-data-simple.csv")
# # subset soildat to 0-10cm depth
# texture <- texture[which(texture$Depth=="0-10 cm"),] 
# soildat$sand[1:80] <- c(rep(texture$sand, each=5))
# soildat$silt[1:80] <- c(rep(texture$silt, each=5))
# soildat$clay[1:80] <- c(rep(texture$clay, each=5))
# soildat$POC[1:80] <- c(rep(texture$POC_baseline, each=5))
# soildat$P[1:80] <- c(rep(texture$P_baseline, each=5))
# soildat$K[1:80] <- c(rep(texture$K_baseline, each=5))
# soildat$Ca[1:80] <- c(rep(texture$Ca_baseline, each=5))
# soildat$Mg[1:80] <- c(rep(texture$Mg_baseline, each=5))
# soildat$Mn[1:80] <- c(rep(texture$Mn_baseline, each=5))
# soildat$Zn[1:80] <- c(rep(texture$Zn_baseline, each=5))
# soildat$Al[1:80] <- c(rep(texture$Al_baseline, each=5))
# soildat$B[1:80] <- c(rep(texture$B_baseline, each=5))
# soildat$Cu[1:80] <- c(rep(texture$Cu_baseline, each=5))
# soildat$Fe[1:80] <- c(rep(texture$Fe_baseline, each=5))
# soildat$Na[1:80] <- c(rep(texture$Na_baseline, each=5))
# soildat$Ni[1:80] <- c(rep(texture$Ni_baseline, each=5))
# soildat$Si[1:80] <- c(rep(texture$Si_baseline, each=5))





################  Bacteria ################  

## there are 13 sampling points:
# fall 2020, (baseline)
# spring 2021, summer 2021, fall 2021, 
# spring 2022, summer 2022, fall 2022
# spring 2023, summer 2023, fall 2023
# spring 2024, summer 2024, fall 2024

sampdate <- c("All_agmicrobiome") 
region <- "16S"
name <- "Bacteria"
# number of sampling periods for which we have sequence data
nsamps <- 13







####################################
### Visualize communities: NMDS
####################################


# select colors for NMDS
cropsyscols <- c("darkgoldenrod4", "darkolivegreen", "darkolivegreen3", "deepskyblue3")
covercols <- c("chocolate", "goldenrod", "forestgreen", "darkcyan", "darkorchid4")
ssnyrcols <- c("white", 
               RColorBrewer::brewer.pal(n = 4, name = "Purples")[2:4], 
               RColorBrewer::brewer.pal(n = 4, name = "YlGn")[2:4], 
               RColorBrewer::brewer.pal(n = 4, name = "Reds")[2:4], 
               RColorBrewer::brewer.pal(n = 4, name = "Blues")[2:4])
yrcols <- c("white", "purple", "aquamarine3", "brown2", "blue")
ssncols <- c("purple", "aquamarine3", "blue")



# import data for abundances and taxonomy
taxabund <- read.csv(paste0("Raw-data/sequence/",sampdate,"/",name,"/seqtab/6seqtab.csv"))
# import data for diversity
divdat <- read.csv(paste0("Processed-data/sequences/",region,"-diversity.csv"))
# import data for diversity
divmat_pa <- read.csv(paste0("Raw-data/sequence/",sampdate,"/divmat_pa_", region,".csv"), row.names = 1)


divdat$Season.Year <- as.factor(divdat$Season.Year)
divdat$Season.Year <- factor(divdat$Season.Year, levels(divdat$Season.Year)[c(1,6,10,2,7,11,3,8,12,4,9,13,5)])
divdat$Season <- as.factor(divdat$Season )
divdat$Season <- factor(divdat$Season, levels(divdat$Season)[c(2,3, 1)])
divdat$Cover <- as.factor(divdat$Cover)
divdat$Cover <- factor(divdat$Cover, levels(divdat$Cover)[c(2,1,4,5,3)])
divdat$Cropping.system <- as.factor(divdat$Cropping.system)
divdat$Cropping.system <- factor(divdat$Cropping.system, levels(divdat$Cropping.system)[c(1,4,3,2)])
divdat$Replicate <- as.factor(divdat$Replicate)
divdat$Year <- as.factor(divdat$Year)
dat_plot <- divdat[which(is.na(divdat$richness)==FALSE),]





##### Community visualization: NMDS of all data

#rem <- c("318_Sum_2021", "220_Fall_2021") #, "210_Spr_2022"
m<-metaMDS(divmat_pa, display="sites", wascores="TRUE", maxit=5000) # [-which(row.names(divmat_pa) %in% rem),]
m
#Stress:     0.2002651 
s<-scores(m, display = c("sites", "species"), shrink = FALSE)
#s
# examine sites... there are some weird ones
# View(s[["sites"]])
# rem <- c("318_Sum_2021", "220_Fall_2021")
distance<-vegdist(divmat_pa, method="bray") # [-which(row.names(divmat_pa) %in% rem),]
ax1_dist<-vegdist(m$points[,1],method="euclidean")
c1 <- cor(distance,ax1_dist)^2 # 
ax2_dist<-vegdist(m$points[,2],method="euclidean")
c2 <- cor(distance,ax2_dist)^2 # 
with(divdat, levels(Season.Year))
scl <-2
#View(s$sites)

xlims <- c(-2.5,3)
ylims <- c(-2,1.5)

# by year
png(paste0("Figures/microbial-diversity/nmds/all-data/",region,"_yr.png"), width=800, height=850, res = 150)
plot(s$sites, type="n", xlim=xlims, ylim=ylims, xlab=bquote(NMDS1~"("*.(round(c1*100,0))*"%)"), ylab=bquote(NMDS2~"("*.(round(c2*100,0))*"%)"), las=1) #  , ylim=c(-1.1,1.5)
with(dat_plot, points(s$sites,   pch =21, bg = alpha(yrcols[Year], 0.9))) # [-which(divdat$new.name %in% rem),]
with(dat_plot, legend(x="topleft", legend=levels(Year), bty="n", col="black", pch=21, pt.cex=1.5, pt.bg=alpha(yrcols,0.8))) # [-which(divdat$new.name %in% rem),]
ordiellipse(m, group=dat_plot$Year, col = alpha(yrcols, 0.9), label=F, draw="polygon", border = "black" ) # [-which(divdat$new.name %in% rem)]
dev.off()

# by season/year
png(paste0("Figures/microbial-diversity/nmds/all-data/",region,"_ssnyr.png"), width=1000, height=1000, res = 150)
plot(s$sites, type="n", xlim=xlims, ylim=ylims, xlab=bquote(NMDS1~"("*.(round(c1*100,0))*"%)"), ylab=bquote(NMDS2~"("*.(round(c2*100,0))*"%)"), las=1) #  , ylim=c(-1.1,1.5)
with(dat_plot, points(s$sites,  scaling =scl, pch = 21, bg = alpha(ssnyrcols[Season.Year], 0.9))) # [-which(divdat$new.name %in% rem),]
with(dat_plot, legend("topleft", legend=levels(Season.Year), bty="n", col="black", pch=21, pt.cex=1.5, pt.bg=alpha(ssnyrcols,0.8))) # [-which(divdat$new.name %in% rem),]
#with(dat_plot, legend("bottomleft", legend=levels(Cropping.system), bty="n", col="black", pch=c(21:24), pt.cex=1.5, pt.bg="white")) # [-which(divdat$new.name %in% rem),]
ordiellipse(m, group=dat_plot$Season.Year, col = alpha(ssnyrcols, 0.9), label=F, draw="polygon", border = "black"  ) # [-which(divdat$new.name %in% rem)]
dev.off()

# by season
png(paste0("Figures/microbial-diversity/nmds/all-data/",region,"_ssn.png"), width=1000, height=1000, res = 150)
plot(s$sites, type="n", xlim=xlims, ylim=ylims, xlab=bquote(NMDS1~"("*.(round(c1*100,0))*"%)"), ylab=bquote(NMDS2~"("*.(round(c2*100,0))*"%)"), las=1) #  , ylim=c(-1.1,1.5)
with(dat_plot, points(s$sites,  scaling =scl, pch = 21, bg = alpha(ssncols[Season], 0.9))) # [-which(divdat$new.name %in% rem),]
with(dat_plot, legend("topleft", legend=levels(Season), bty="n", col="black", pch=21, pt.cex=1.5, pt.bg=alpha(ssncols,0.8))) # [-which(divdat$new.name %in% rem),]
#with(dat_plot, legend("bottomleft", legend=levels(Cropping.system), bty="n", col="black", pch=c(21:24), pt.cex=1.5, pt.bg="white")) # [-which(divdat$new.name %in% rem),]
ordiellipse(m, group=dat_plot$Season, col = alpha(ssncols, 0.9), label=F, draw="polygon", border = "black"  ) # [-which(divdat$new.name %in% rem)]
dev.off()

# by cover
png(paste0("Figures/microbial-diversity/nmds/all-data/",region,"_cover.png"), width=600, height=650, res = 150)
plot(s$sites, type="n", xlim=xlims, ylim=ylims,  xlab=bquote(NMDS1~"("*.(round(c1*100,0))*"%)"), ylab=bquote(NMDS2~"("*.(round(c2*100,0))*"%)"), las=1)
with(dat_plot, points(s$sites,  scaling =scl, pch =21, bg = alpha(covercols[Cover], 0.9)))
with(dat_plot, legend("bottomleft", legend=levels(Cover), bty="n", col="black", pch=21, pt.cex=1.5, pt.bg=alpha(covercols,0.8)))
ordiellipse(m, group=dat_plot$Cover, col = alpha(covercols, 0.9), label=F, draw="polygon", border = "black"  )
dev.off()

# by cropping system
png(paste0("Figures/microbial-diversity/nmds/all-data/",region,"_cropsys.png"), width=700, height=750, res = 150)
plot(s$sites, type="n", xlim=xlims, ylim=ylims, xlab=bquote(NMDS1~"("*.(round(c1*100,0))*"%)"), ylab=bquote(NMDS2~"("*.(round(c2*100,0))*"%)"), las=1)
with(dat_plot, points(s$sites,  scaling =scl, pch =21, bg = alpha(cropsyscols[Cropping.system], 0.9)))
with(dat_plot, legend("bottomleft", legend=levels(Cropping.system), bty="n", col="black", pch=21, pt.cex=1.5, pt.bg=alpha(cropsyscols,0.8)))
ordiellipse(m, group=dat_plot$Cropping.system, col = alpha(cropsyscols, 0.9), label=F, draw="polygon", border = "black"  )
dev.off()





##### NMDS: subset to each season

# Spring 2021
ssn <- "Spring 2021"
xlims <- c(-0.8,0.8)
divdat_2 <- dat_plot[which(dat_plot$Season.Year==ssn),]
divmat_pa_2 <- divmat_pa[which(dat_plot$Season.Year==ssn),]
#rem <- c("318_Sum_2021", "220_Fall_2021")
m<-metaMDS(divmat_pa_2, display="sites", wascores="TRUE", maxit=5000)
s<-scores(m, display = c("sites", "species"), shrink = FALSE)
distance<-vegdist(divmat_pa_2, method="bray")
ax1_dist<-vegdist(m$points[,1],method="euclidean")
c1<-cor(distance,ax1_dist)^2 # 0.3852263
ax2_dist<-vegdist(m$points[,2],method="euclidean")
c2<-cor(distance,ax2_dist)^2 # 0.3230541
with(divdat_2, levels(Season.Year))
scl <-2
# by cover
png(paste0("Figures/microbial-diversity/nmds/by-season/",region, "_", ssn, "_cover.png"),  width=700, height=750, res = 150)
plot(s$sites, xlim=xlims, type="n", main=bquote(.(ssn)), xlab=bquote(NMDS1~"("*.(round(c1*100,0))*"%)"), ylab=bquote(NMDS2~"("*.(round(c2*100,0))*"%)"), las=1)#, xlim=c(-1,2.5), ylim=c(-1.1,1.5))
with(divdat_2, points(s$sites,  scaling =scl, pch =21, bg = alpha(covercols[Cover], 0.9)))
with(divdat_2, legend("topleft", legend=levels(Cover), bty="n", col="black", pch=21, pt.cex=1.5, pt.bg=alpha(covercols,0.9)))
ordiellipse(m, group=divdat_2$Cover, col = alpha(covercols, 0.9), label=F, draw="polygon", border = "black" )
dev.off()
# by cropping system
png(paste0("Figures/microbial-diversity/nmds/by-season/",region,"_", ssn, "_cropsys.png"), width=700, height=750, res = 150)
plot(s$sites, xlim=xlims, type="n",main=bquote(.(ssn)), xlab=bquote(NMDS1~"("*.(round(c1*100,0))*"%)"), ylab=bquote(NMDS2~"("*.(round(c2*100,0))*"%)"), las=1)
with(divdat_2, points(s$sites,  scaling =scl, pch =21, bg = alpha(cropsyscols[Cropping.system], 0.9)))
with(divdat_2, legend("topleft", legend=levels(Cropping.system), bty="n", col="black", pch=21, pt.cex=1.5, pt.bg=alpha(cropsyscols,0.9)))
ordiellipse(m, group=divdat_2$Cropping.system, col = alpha(cropsyscols, 0.9), label=F, draw="polygon", border = "black" )
dev.off()

# Summer 2021
ssn <- "Summer 2021"
#xlims <- c(-0.8,0.4)
divdat_2 <- dat_plot[which(dat_plot$Season.Year==ssn),]
divmat_pa_2 <- divmat_pa[which(dat_plot$Season.Year==ssn),]
#rem <- c("318_Sum_2021", "220_Fall_2021")
m<-metaMDS(divmat_pa_2, display="sites", wascores="TRUE", maxit=5000)
s<-scores(m, display = c("sites", "species"), shrink = FALSE)
distance<-vegdist(divmat_pa_2, method="bray")
ax1_dist<-vegdist(m$points[,1],method="euclidean")
c1<-cor(distance,ax1_dist)^2 # 0.3852263
ax2_dist<-vegdist(m$points[,2],method="euclidean")
c2<-cor(distance,ax2_dist)^2 # 0.3230541
with(divdat_2, levels(Season.Year))
scl <-2
# by cover
png(paste0("Figures/microbial-diversity/nmds/by-season/",region, "_", ssn, "_cover.png"),  width=700, height=750, res = 150)
plot(s$sites, type="n", main=bquote(.(ssn)), xlab=bquote(NMDS1~"("*.(round(c1*100,0))*"%)"), ylab=bquote(NMDS2~"("*.(round(c2*100,0))*"%)"), las=1)#, xlim=c(-1,2.5), ylim=c(-1.1,1.5))
with(divdat_2, points(s$sites,  scaling =scl, pch =21, bg = alpha(covercols[Cover], 0.9)))
#with(divdat_2, legend("topleft", legend=levels(Cover), bty="n", col="black", pch=21, pt.cex=1.5, pt.bg=alpha(covercols,0.9)))
ordiellipse(m, group=divdat_2$Cover, col = alpha(covercols, 0.9), label=F, draw="polygon", border = "black" )
dev.off()
# by cropping system
png(paste0("Figures/microbial-diversity/nmds/by-season/",region,"_", ssn, "_cropsys.png"), width=700, height=750, res = 150)
plot(s$sites, type="n",main=bquote(.(ssn)), xlab=bquote(NMDS1~"("*.(round(c1*100,0))*"%)"), ylab=bquote(NMDS2~"("*.(round(c2*100,0))*"%)"), las=1)
with(divdat_2, points(s$sites,  scaling =scl, pch =21, bg = alpha(cropsyscols[Cropping.system], 0.9)))
#with(divdat_2, legend("topleft", legend=levels(Cropping.system), bty="n", col="black", pch=21, pt.cex=1.5, pt.bg=alpha(cropsyscols,0.9)))
ordiellipse(m, group=divdat_2$Cropping.system, col = alpha(cropsyscols, 0.9), label=F, draw="polygon", border = "black" )
dev.off()

# Fall 2021
ssn <- "Fall 2021"
#xlims <- c(-0.8,0.4)
divdat_2 <- dat_plot[which(dat_plot$Season.Year==ssn),]
divmat_pa_2 <- divmat_pa[which(dat_plot$Season.Year==ssn),]
#rem <- c("318_Sum_2021", "220_Fall_2021")
m<-metaMDS(divmat_pa_2, display="sites", wascores="TRUE", maxit=5000)
s<-scores(m, display = c("sites", "species"), shrink = FALSE)
distance<-vegdist(divmat_pa_2, method="bray")
ax1_dist<-vegdist(m$points[,1],method="euclidean")
c1<-cor(distance,ax1_dist)^2 # 0.3852263
ax2_dist<-vegdist(m$points[,2],method="euclidean")
c2<-cor(distance,ax2_dist)^2 # 0.3230541
with(divdat_2, levels(Season.Year))
scl <-2
# by cover
png(paste0("Figures/microbial-diversity/nmds/by-season/",region, "_", ssn, "_cover.png"),  width=700, height=750, res = 150)
plot(s$sites, type="n", main=bquote(.(ssn)), xlab=bquote(NMDS1~"("*.(round(c1*100,0))*"%)"), ylab=bquote(NMDS2~"("*.(round(c2*100,0))*"%)"), las=1)#, xlim=c(-1,2.5), ylim=c(-1.1,1.5))
with(divdat_2, points(s$sites,  scaling =scl, pch =21, bg = alpha(covercols[Cover], 0.9)))
#with(divdat_2, legend("topleft", legend=levels(Cover), bty="n", col="black", pch=21, pt.cex=1.5, pt.bg=alpha(covercols,0.9)))
ordiellipse(m, group=divdat_2$Cover, col = alpha(covercols, 0.9), label=F, draw="polygon", border = "black" )
dev.off()
# by cropping system
png(paste0("Figures/microbial-diversity/nmds/by-season/",region,"_", ssn, "_cropsys.png"), width=700, height=750, res = 150)
plot(s$sites, type="n",main=bquote(.(ssn)), xlab=bquote(NMDS1~"("*.(round(c1*100,0))*"%)"), ylab=bquote(NMDS2~"("*.(round(c2*100,0))*"%)"), las=1)
with(divdat_2, points(s$sites,  scaling =scl, pch =21, bg = alpha(cropsyscols[Cropping.system], 0.9)))
#with(divdat_2, legend("topleft", legend=levels(Cropping.system), bty="n", col="black", pch=21, pt.cex=1.5, pt.bg=alpha(cropsyscols,0.9)))
ordiellipse(m, group=divdat_2$Cropping.system, col = alpha(cropsyscols, 0.9), label=F, draw="polygon", border = "black" )
dev.off()





# Spring 2022
ssn <- "Spring 2022"
#xlims <- c(-0.8,0.8)
divdat_2 <- dat_plot[which(dat_plot$Season.Year==ssn),]
divmat_pa_2 <- divmat_pa[which(dat_plot$Season.Year==ssn),]
#rem <- c("318_Sum_2022", "220_Fall_2022")
m<-metaMDS(divmat_pa_2, display="sites", wascores="TRUE", maxit=5000)
s<-scores(m, display = c("sites", "species"), shrink = FALSE)
distance<-vegdist(divmat_pa_2, method="bray")
ax1_dist<-vegdist(m$points[,1],method="euclidean")
c1<-cor(distance,ax1_dist)^2 # 0.3852263
ax2_dist<-vegdist(m$points[,2],method="euclidean")
c2<-cor(distance,ax2_dist)^2 # 0.3230541
with(divdat_2, levels(Season.Year))
scl <-2
# by cover
png(paste0("Figures/microbial-diversity/nmds/by-season/",region, "_", ssn, "_cover.png"),  width=700, height=750, res = 150)
plot(s$sites,  type="n", main=bquote(.(ssn)), xlab=bquote(NMDS1~"("*.(round(c1*100,0))*"%)"), ylab=bquote(NMDS2~"("*.(round(c2*100,0))*"%)"), las=1)#, xlim=c(-1,2.5), ylim=c(-1.1,1.5))
with(divdat_2, points(s$sites,  scaling =scl, pch =21, bg = alpha(covercols[Cover], 0.9)))
#with(divdat_2, legend("topleft", legend=levels(Cover), bty="n", col="black", pch=21, pt.cex=1.5, pt.bg=alpha(covercols,0.9)))
ordiellipse(m, group=divdat_2$Cover, col = alpha(covercols, 0.9), label=F, draw="polygon", border = "black" )
dev.off()
# by cropping system
png(paste0("Figures/microbial-diversity/nmds/by-season/",region,"_", ssn, "_cropsys.png"), width=700, height=750, res = 150)
plot(s$sites, type="n",main=bquote(.(ssn)), xlab=bquote(NMDS1~"("*.(round(c1*100,0))*"%)"), ylab=bquote(NMDS2~"("*.(round(c2*100,0))*"%)"), las=1)
with(divdat_2, points(s$sites,  scaling =scl, pch =21, bg = alpha(cropsyscols[Cropping.system], 0.9)))
#with(divdat_2, legend("topleft", legend=levels(Cropping.system), bty="n", col="black", pch=21, pt.cex=1.5, pt.bg=alpha(cropsyscols,0.9)))
ordiellipse(m, group=divdat_2$Cropping.system, col = alpha(cropsyscols, 0.9), label=F, draw="polygon", border = "black" )
dev.off()

# Summer 2022
ssn <- "Summer 2022"
#xlims <- c(-0.8,0.4)
divdat_2 <- dat_plot[which(dat_plot$Season.Year==ssn),]
divmat_pa_2 <- divmat_pa[which(dat_plot$Season.Year==ssn),]
#rem <- c("318_Sum_2022", "220_Fall_2022")
m<-metaMDS(divmat_pa_2, display="sites", wascores="TRUE", maxit=5000)
s<-scores(m, display = c("sites", "species"), shrink = FALSE)
distance<-vegdist(divmat_pa_2, method="bray")
ax1_dist<-vegdist(m$points[,1],method="euclidean")
c1<-cor(distance,ax1_dist)^2 # 0.3852263
ax2_dist<-vegdist(m$points[,2],method="euclidean")
c2<-cor(distance,ax2_dist)^2 # 0.3230541
with(divdat_2, levels(Season.Year))
scl <-2
# by cover
png(paste0("Figures/microbial-diversity/nmds/by-season/",region, "_", ssn, "_cover.png"),  width=700, height=750, res = 150)
plot(s$sites, type="n", main=bquote(.(ssn)), xlab=bquote(NMDS1~"("*.(round(c1*100,0))*"%)"), ylab=bquote(NMDS2~"("*.(round(c2*100,0))*"%)"), las=1)#, xlim=c(-1,2.5), ylim=c(-1.1,1.5))
with(divdat_2, points(s$sites,  scaling =scl, pch =21, bg = alpha(covercols[Cover], 0.9)))
#with(divdat_2, legend("topleft", legend=levels(Cover), bty="n", col="black", pch=21, pt.cex=1.5, pt.bg=alpha(covercols,0.9)))
ordiellipse(m, group=divdat_2$Cover, col = alpha(covercols, 0.9), label=F, draw="polygon", border = "black" )
dev.off()
# by cropping system
png(paste0("Figures/microbial-diversity/nmds/by-season/",region,"_", ssn, "_cropsys.png"), width=700, height=750, res = 150)
plot(s$sites, type="n",main=bquote(.(ssn)), xlab=bquote(NMDS1~"("*.(round(c1*100,0))*"%)"), ylab=bquote(NMDS2~"("*.(round(c2*100,0))*"%)"), las=1)
with(divdat_2, points(s$sites,  scaling =scl, pch =21, bg = alpha(cropsyscols[Cropping.system], 0.9)))
#with(divdat_2, legend("topleft", legend=levels(Cropping.system), bty="n", col="black", pch=21, pt.cex=1.5, pt.bg=alpha(cropsyscols,0.9)))
ordiellipse(m, group=divdat_2$Cropping.system, col = alpha(cropsyscols, 0.9), label=F, draw="polygon", border = "black" )
dev.off()

# Fall 2022
ssn <- "Fall 2022"
#xlims <- c(-0.8,0.4)
divdat_2 <- dat_plot[which(dat_plot$Season.Year==ssn),]
divmat_pa_2 <- divmat_pa[which(dat_plot$Season.Year==ssn),]
#rem <- c("318_Sum_2022", "220_Fall_2022")
m<-metaMDS(divmat_pa_2, display="sites", wascores="TRUE", maxit=5000)
s<-scores(m, display = c("sites", "species"), shrink = FALSE)
distance<-vegdist(divmat_pa_2, method="bray")
ax1_dist<-vegdist(m$points[,1],method="euclidean")
c1<-cor(distance,ax1_dist)^2 # 0.3852263
ax2_dist<-vegdist(m$points[,2],method="euclidean")
c2<-cor(distance,ax2_dist)^2 # 0.3230541
with(divdat_2, levels(Season.Year))
scl <-2
# by cover
png(paste0("Figures/microbial-diversity/nmds/by-season/",region, "_", ssn, "_cover.png"),  width=700, height=750, res = 150)
plot(s$sites, type="n", main=bquote(.(ssn)), xlab=bquote(NMDS1~"("*.(round(c1*100,0))*"%)"), ylab=bquote(NMDS2~"("*.(round(c2*100,0))*"%)"), las=1)#, xlim=c(-1,2.5), ylim=c(-1.1,1.5))
with(divdat_2, points(s$sites,  scaling =scl, pch =21, bg = alpha(covercols[Cover], 0.9)))
#with(divdat_2, legend("topleft", legend=levels(Cover), bty="n", col="black", pch=21, pt.cex=1.5, pt.bg=alpha(covercols,0.9)))
ordiellipse(m, group=divdat_2$Cover, col = alpha(covercols, 0.9), label=F, draw="polygon", border = "black" )
dev.off()
# by cropping system
png(paste0("Figures/microbial-diversity/nmds/by-season/",region,"_", ssn, "_cropsys.png"), width=700, height=750, res = 150)
plot(s$sites, type="n",main=bquote(.(ssn)), xlab=bquote(NMDS1~"("*.(round(c1*100,0))*"%)"), ylab=bquote(NMDS2~"("*.(round(c2*100,0))*"%)"), las=1)
with(divdat_2, points(s$sites,  scaling =scl, pch =21, bg = alpha(cropsyscols[Cropping.system], 0.9)))
#with(divdat_2, legend("topleft", legend=levels(Cropping.system), bty="n", col="black", pch=21, pt.cex=1.5, pt.bg=alpha(cropsyscols,0.9)))
ordiellipse(m, group=divdat_2$Cropping.system, col = alpha(cropsyscols, 0.9), label=F, draw="polygon", border = "black" )
dev.off()






# Spring 2023
ssn <- "Spring 2023"
#xlims <- c(-0.8,0.8)
divdat_2 <- dat_plot[which(dat_plot$Season.Year==ssn),]
divmat_pa_2 <- divmat_pa[which(dat_plot$Season.Year==ssn),]
#rem <- c("318_Sum_2023", "220_Fall_2023")
m<-metaMDS(divmat_pa_2, display="sites", wascores="TRUE", maxit=5000)
s<-scores(m, display = c("sites", "species"), shrink = FALSE)
distance<-vegdist(divmat_pa_2, method="bray")
ax1_dist<-vegdist(m$points[,1],method="euclidean")
c1<-cor(distance,ax1_dist)^2 # 0.3852263
ax2_dist<-vegdist(m$points[,2],method="euclidean")
c2<-cor(distance,ax2_dist)^2 # 0.3230541
with(divdat_2, levels(Season.Year))
scl <-2
# by cover
png(paste0("Figures/microbial-diversity/nmds/by-season/",region, "_", ssn, "_cover.png"),  width=700, height=750, res = 150)
plot(s$sites,  type="n", main=bquote(.(ssn)), xlab=bquote(NMDS1~"("*.(round(c1*100,0))*"%)"), ylab=bquote(NMDS2~"("*.(round(c2*100,0))*"%)"), las=1)#, xlim=c(-1,2.5), ylim=c(-1.1,1.5))
with(divdat_2, points(s$sites,  scaling =scl, pch =21, bg = alpha(covercols[Cover], 0.9)))
#with(divdat_2, legend("topleft", legend=levels(Cover), bty="n", col="black", pch=21, pt.cex=1.5, pt.bg=alpha(covercols,0.9)))
ordiellipse(m, group=divdat_2$Cover, col = alpha(covercols, 0.9), label=F, draw="polygon", border = "black" )
dev.off()
# by cropping system
png(paste0("Figures/microbial-diversity/nmds/by-season/",region,"_", ssn, "_cropsys.png"), width=700, height=750, res = 150)
plot(s$sites, type="n",main=bquote(.(ssn)), xlab=bquote(NMDS1~"("*.(round(c1*100,0))*"%)"), ylab=bquote(NMDS2~"("*.(round(c2*100,0))*"%)"), las=1)
with(divdat_2, points(s$sites,  scaling =scl, pch =21, bg = alpha(cropsyscols[Cropping.system], 0.9)))
#with(divdat_2, legend("topleft", legend=levels(Cropping.system), bty="n", col="black", pch=21, pt.cex=1.5, pt.bg=alpha(cropsyscols,0.9)))
ordiellipse(m, group=divdat_2$Cropping.system, col = alpha(cropsyscols, 0.9), label=F, draw="polygon", border = "black" )
dev.off()

# Summer 2023
ssn <- "Summer 2023"
#xlims <- c(-0.8,0.4)
divdat_2 <- dat_plot[which(dat_plot$Season.Year==ssn),]
divmat_pa_2 <- divmat_pa[which(dat_plot$Season.Year==ssn),]
#rem <- c("318_Sum_2023", "220_Fall_2023")
m<-metaMDS(divmat_pa_2, display="sites", wascores="TRUE", maxit=5000)
s<-scores(m, display = c("sites", "species"), shrink = FALSE)
distance<-vegdist(divmat_pa_2, method="bray")
ax1_dist<-vegdist(m$points[,1],method="euclidean")
c1<-cor(distance,ax1_dist)^2 # 0.3852263
ax2_dist<-vegdist(m$points[,2],method="euclidean")
c2<-cor(distance,ax2_dist)^2 # 0.3230541
with(divdat_2, levels(Season.Year))
scl <-2
# by cover
png(paste0("Figures/microbial-diversity/nmds/by-season/",region, "_", ssn, "_cover.png"),  width=700, height=750, res = 150)
plot(s$sites, type="n", main=bquote(.(ssn)), xlab=bquote(NMDS1~"("*.(round(c1*100,0))*"%)"), ylab=bquote(NMDS2~"("*.(round(c2*100,0))*"%)"), las=1)#, xlim=c(-1,2.5), ylim=c(-1.1,1.5))
with(divdat_2, points(s$sites,  scaling =scl, pch =21, bg = alpha(covercols[Cover], 0.9)))
#with(divdat_2, legend("topleft", legend=levels(Cover), bty="n", col="black", pch=21, pt.cex=1.5, pt.bg=alpha(covercols,0.9)))
ordiellipse(m, group=divdat_2$Cover, col = alpha(covercols, 0.9), label=F, draw="polygon", border = "black" )
dev.off()
# by cropping system
png(paste0("Figures/microbial-diversity/nmds/by-season/",region,"_", ssn, "_cropsys.png"), width=700, height=750, res = 150)
plot(s$sites, type="n",main=bquote(.(ssn)), xlab=bquote(NMDS1~"("*.(round(c1*100,0))*"%)"), ylab=bquote(NMDS2~"("*.(round(c2*100,0))*"%)"), las=1)
with(divdat_2, points(s$sites,  scaling =scl, pch =21, bg = alpha(cropsyscols[Cropping.system], 0.9)))
#with(divdat_2, legend("topleft", legend=levels(Cropping.system), bty="n", col="black", pch=21, pt.cex=1.5, pt.bg=alpha(cropsyscols,0.9)))
ordiellipse(m, group=divdat_2$Cropping.system, col = alpha(cropsyscols, 0.9), label=F, draw="polygon", border = "black" )
dev.off()

# Fall 2023
ssn <- "Fall 2023"
#xlims <- c(-0.8,0.4)
divdat_2 <- dat_plot[which(dat_plot$Season.Year==ssn),]
divmat_pa_2 <- divmat_pa[which(dat_plot$Season.Year==ssn),]
#rem <- c("318_Sum_2023", "220_Fall_2023")
m<-metaMDS(divmat_pa_2, display="sites", wascores="TRUE", maxit=5000)
s<-scores(m, display = c("sites", "species"), shrink = FALSE)
distance<-vegdist(divmat_pa_2, method="bray")
ax1_dist<-vegdist(m$points[,1],method="euclidean")
c1<-cor(distance,ax1_dist)^2 # 0.3852263
ax2_dist<-vegdist(m$points[,2],method="euclidean")
c2<-cor(distance,ax2_dist)^2 # 0.3230541
with(divdat_2, levels(Season.Year))
scl <-2
# by cover
png(paste0("Figures/microbial-diversity/nmds/by-season/",region, "_", ssn, "_cover.png"),  width=700, height=750, res = 150)
plot(s$sites, type="n", main=bquote(.(ssn)), xlab=bquote(NMDS1~"("*.(round(c1*100,0))*"%)"), ylab=bquote(NMDS2~"("*.(round(c2*100,0))*"%)"), las=1)#, xlim=c(-1,2.5), ylim=c(-1.1,1.5))
with(divdat_2, points(s$sites,  scaling =scl, pch =21, bg = alpha(covercols[Cover], 0.9)))
#with(divdat_2, legend("topleft", legend=levels(Cover), bty="n", col="black", pch=21, pt.cex=1.5, pt.bg=alpha(covercols,0.9)))
ordiellipse(m, group=divdat_2$Cover, col = alpha(covercols, 0.9), label=F, draw="polygon", border = "black" )
dev.off()
# by cropping system
png(paste0("Figures/microbial-diversity/nmds/by-season/",region,"_", ssn, "_cropsys.png"), width=700, height=750, res = 150)
plot(s$sites, type="n",main=bquote(.(ssn)), xlab=bquote(NMDS1~"("*.(round(c1*100,0))*"%)"), ylab=bquote(NMDS2~"("*.(round(c2*100,0))*"%)"), las=1)
with(divdat_2, points(s$sites,  scaling =scl, pch =21, bg = alpha(cropsyscols[Cropping.system], 0.9)))
#with(divdat_2, legend("topleft", legend=levels(Cropping.system), bty="n", col="black", pch=21, pt.cex=1.5, pt.bg=alpha(cropsyscols,0.9)))
ordiellipse(m, group=divdat_2$Cropping.system, col = alpha(cropsyscols, 0.9), label=F, draw="polygon", border = "black" )
dev.off()






# Spring 2024
ssn <- "Spring 2024"
#xlims <- c(-0.8,0.8)
divdat_2 <- dat_plot[which(dat_plot$Season.Year==ssn),]
divmat_pa_2 <- divmat_pa[which(dat_plot$Season.Year==ssn),]
#rem <- c("318_Sum_2024", "220_Fall_2024")
m<-metaMDS(divmat_pa_2, display="sites", wascores="TRUE", maxit=5000)
s<-scores(m, display = c("sites", "species"), shrink = FALSE)
distance<-vegdist(divmat_pa_2, method="bray")
ax1_dist<-vegdist(m$points[,1],method="euclidean")
c1<-cor(distance,ax1_dist)^2 # 0.3852263
ax2_dist<-vegdist(m$points[,2],method="euclidean")
c2<-cor(distance,ax2_dist)^2 # 0.3230541
with(divdat_2, levels(Season.Year))
scl <-2
# by cover
png(paste0("Figures/microbial-diversity/nmds/by-season/",region, "_", ssn, "_cover.png"),  width=700, height=750, res = 150)
plot(s$sites,  type="n", main=bquote(.(ssn)), xlab=bquote(NMDS1~"("*.(round(c1*100,0))*"%)"), ylab=bquote(NMDS2~"("*.(round(c2*100,0))*"%)"), las=1)#, xlim=c(-1,2.5), ylim=c(-1.1,1.5))
with(divdat_2, points(s$sites,  scaling =scl, pch =21, bg = alpha(covercols[Cover], 0.9)))
#with(divdat_2, legend("topleft", legend=levels(Cover), bty="n", col="black", pch=21, pt.cex=1.5, pt.bg=alpha(covercols,0.9)))
ordiellipse(m, group=divdat_2$Cover, col = alpha(covercols, 0.9), label=F, draw="polygon", border = "black" )
dev.off()
# by cropping system
png(paste0("Figures/microbial-diversity/nmds/by-season/",region,"_", ssn, "_cropsys.png"), width=700, height=750, res = 150)
plot(s$sites, type="n",main=bquote(.(ssn)), xlab=bquote(NMDS1~"("*.(round(c1*100,0))*"%)"), ylab=bquote(NMDS2~"("*.(round(c2*100,0))*"%)"), las=1)
with(divdat_2, points(s$sites,  scaling =scl, pch =21, bg = alpha(cropsyscols[Cropping.system], 0.9)))
#with(divdat_2, legend("topleft", legend=levels(Cropping.system), bty="n", col="black", pch=21, pt.cex=1.5, pt.bg=alpha(cropsyscols,0.9)))
ordiellipse(m, group=divdat_2$Cropping.system, col = alpha(cropsyscols, 0.9), label=F, draw="polygon", border = "black" )
dev.off()

# Summer 2024 (remove 203_Summer_2024)
ssn <- "Summer 2024"
#xlims <- c(-0.8,0.4)
divdat_2 <- dat_plot[which(dat_plot$Season.Year==ssn),]
divmat_pa_2 <- divmat_pa[which(dat_plot$Season.Year==ssn),]
#rem <- c("318_Sum_2024", "220_Fall_2024")
m<-metaMDS(divmat_pa_2, display="sites", wascores="TRUE", maxit=5000)
s<-scores(m, display = c("sites", "species"), shrink = FALSE)
distance<-vegdist(divmat_pa_2, method="bray")
ax1_dist<-vegdist(m$points[,1],method="euclidean")
c1<-cor(distance,ax1_dist)^2 # 0.3852263
ax2_dist<-vegdist(m$points[,2],method="euclidean")
c2<-cor(distance,ax2_dist)^2 # 0.3230541
with(divdat_2, levels(Season.Year))
scl <-2
# by cover
png(paste0("Figures/microbial-diversity/nmds/by-season/",region, "_", ssn, "_cover.png"),  width=700, height=750, res = 150)
plot(s$sites, type="n", main=bquote(.(ssn)), xlab=bquote(NMDS1~"("*.(round(c1*100,0))*"%)"), ylab=bquote(NMDS2~"("*.(round(c2*100,0))*"%)"), las=1)#, xlim=c(-1,2.5), ylim=c(-1.1,1.5))
with(divdat_2, points(s$sites,  scaling =scl, pch =21, bg = alpha(covercols[Cover], 0.9)))
#with(divdat_2, legend("topleft", legend=levels(Cover), bty="n", col="black", pch=21, pt.cex=1.5, pt.bg=alpha(covercols,0.9)))
ordiellipse(m, group=divdat_2$Cover, col = alpha(covercols, 0.9), label=F, draw="polygon", border = "black" )
dev.off()
# by cropping system
png(paste0("Figures/microbial-diversity/nmds/by-season/",region,"_", ssn, "_cropsys.png"), width=700, height=750, res = 150)
plot(s$sites, type="n",main=bquote(.(ssn)), xlab=bquote(NMDS1~"("*.(round(c1*100,0))*"%)"), ylab=bquote(NMDS2~"("*.(round(c2*100,0))*"%)"), las=1)
with(divdat_2, points(s$sites,  scaling =scl, pch =21, bg = alpha(cropsyscols[Cropping.system], 0.9)))
#with(divdat_2, legend("topleft", legend=levels(Cropping.system), bty="n", col="black", pch=21, pt.cex=1.5, pt.bg=alpha(cropsyscols,0.9)))
ordiellipse(m, group=divdat_2$Cropping.system, col = alpha(cropsyscols, 0.9), label=F, draw="polygon", border = "black" )
dev.off()

# Fall 2024 (remove 114_Fall_2024)
ssn <- "Fall 2024"
#xlims <- c(-0.8,0.4)
divdat_2 <- dat_plot[which(dat_plot$Season.Year==ssn),]
divmat_pa_2 <- divmat_pa[which(dat_plot$Season.Year==ssn),]
#rem <- c("318_Sum_2024", "220_Fall_2024")
m<-metaMDS(divmat_pa_2, display="sites", wascores="TRUE", maxit=5000)
s<-scores(m, display = c("sites", "species"), shrink = FALSE)
distance<-vegdist(divmat_pa_2, method="bray")
ax1_dist<-vegdist(m$points[,1],method="euclidean")
c1<-cor(distance,ax1_dist)^2 # 0.3852263
ax2_dist<-vegdist(m$points[,2],method="euclidean")
c2<-cor(distance,ax2_dist)^2 # 0.3230541
with(divdat_2, levels(Season.Year))
scl <-2
# by cover
png(paste0("Figures/microbial-diversity/nmds/by-season/",region, "_", ssn, "_cover.png"),  width=700, height=750, res = 150)
plot(s$sites, type="n", main=bquote(.(ssn)), xlab=bquote(NMDS1~"("*.(round(c1*100,0))*"%)"), ylab=bquote(NMDS2~"("*.(round(c2*100,0))*"%)"), las=1)#, xlim=c(-1,2.5), ylim=c(-1.1,1.5))
with(divdat_2, points(s$sites,  scaling =scl, pch =21, bg = alpha(covercols[Cover], 0.9)))
#with(divdat_2, legend("topleft", legend=levels(Cover), bty="n", col="black", pch=21, pt.cex=1.5, pt.bg=alpha(covercols,0.9)))
ordiellipse(m, group=divdat_2$Cover, col = alpha(covercols, 0.9), label=F, draw="polygon", border = "black" )
dev.off()
# by cropping system
png(paste0("Figures/microbial-diversity/nmds/by-season/",region,"_", ssn, "_cropsys.png"), width=700, height=750, res = 150)
plot(s$sites, type="n",main=bquote(.(ssn)), xlab=bquote(NMDS1~"("*.(round(c1*100,0))*"%)"), ylab=bquote(NMDS2~"("*.(round(c2*100,0))*"%)"), las=1)
with(divdat_2, points(s$sites,  scaling =scl, pch =21, bg = alpha(cropsyscols[Cropping.system], 0.9)))
#with(divdat_2, legend("topleft", legend=levels(Cropping.system), bty="n", col="black", pch=21, pt.cex=1.5, pt.bg=alpha(cropsyscols,0.9)))
ordiellipse(m, group=divdat_2$Cropping.system, col = alpha(cropsyscols, 0.9), label=F, draw="polygon", border = "black" )
dev.off()






####################################
### Treatment effects on communities: PERMANOVA
####################################

# import data for abundances and taxonomy
taxabund <- read.csv(paste0("Raw-data/sequence/",sampdate,"/",name,"/seqtab/6seqtab.csv"))
# import data for diversity
divdat <- read.csv(paste0("Processed-data/sequences/", region, "-diversity.csv"))
# import data for diversity
divmat_pa <- read.csv(paste0("Raw-data/sequence/",sampdate,"/divmat_pa_",region,".csv"), row.names = 1)


divdat$Season.Year <- as.factor(divdat$Season.Year)
divdat$Season.Year <- factor(divdat$Season.Year, levels(divdat$Season.Year)[c(1,6,10,2,7,11,3,8,12,4,9,13,5)])
divdat$Season <- as.factor(divdat$Season )
divdat$Season <- factor(divdat$Season, levels(divdat$Season)[c(2,3, 1)])
divdat$Cover <- as.factor(divdat$Cover)
divdat$Cover <- factor(divdat$Cover, levels(divdat$Cover)[c(2,1,4,5,3)])
divdat$Cropping.system <- as.factor(divdat$Cropping.system)
divdat$Cropping.system <- factor(divdat$Cropping.system, levels(divdat$Cropping.system)[c(1,4,3,2)])
divdat$Replicate <- as.factor(divdat$Replicate)
divdat$Year <- as.factor(divdat$Year)


# select data for PERMANOVA
dat_plot <- divdat[which(is.na(divdat$richness)==FALSE),]


### permanova: all data
colnames(divmat_pa)[1:5]; rownames(divmat_pa)[1:5] # asvs as columns, samples as rows
rownames(dat_plot) <- rownames(divmat_pa) # plot number and year
dat_plot$Season.Year 

# effects of treatments: 
permmod <- adonis2(divmat_pa ~ Cover*Cropping.system*Season.Year, 
                   data = dat_plot, 
                   method="bray", 
                   permutations = 199, 
                   by="term", 
                   strata=dat_plot$Replicate)

permmod2 <- round(as.data.frame(permmod), 3)
# ssp <- paste0("F=", permmod2$F, ", Df=", permmod2$Df, ", p=", permmod2$`Pr(>F)`, ", R2=", permmod2$R2)
# mod.effects$result <- ssp
write.csv(permmod2, paste0("Model-output/permanova/",name,"_permanova.csv"))






### permanova: subset to each season
mod.effects <- data.frame(effect=c("Cover","Cropping.system","Cropping.system:Cover","Residual","Total")) 

# Fall 2020
divmat_pa_year2 <- divmat_pa[which(dat_plot$Season.Year=="Fall 2020"),]
divdat_year2 <- dat_plot[which(dat_plot$Season.Year=="Fall 2020"),]
permmod <- adonis2(divmat_pa_year2 ~ Cover*Cropping.system, data = divdat_year2, by = "term",
                   method="bray", permutations = 499, strata = divdat_year2$Replicate)
permmod2 <- round(as.data.frame(permmod), 3)
# append to a mod.effects dataframe
ssp <- paste0("F=", permmod2$F, ", Df=", permmod2$Df, ", p=", permmod2$`Pr(>F)`, ", R2=", permmod2$R2)
mod.effects$'Fall 2020' <- ssp

# Spring 2021
divmat_pa_year2 <- divmat_pa[which(dat_plot$Season.Year=="Spring 2021"),]
dat_plot_year2 <- dat_plot[which(dat_plot$Season.Year=="Spring 2021"),]
permmod <- adonis2(divmat_pa_year2 ~ Cover*Cropping.system, data = dat_plot_year2, by = "term",
                   method="bray", permutations = 499, strata = dat_plot_year2$Replicate)
permmod2 <- round(as.data.frame(permmod), 3)
# append to a mod.effects dataframe
ssp <- paste0("F=", permmod2$F, ", Df=", permmod2$Df, ", p=", permmod2$`Pr(>F)`, ", R2=", permmod2$R2)
mod.effects$'Spring 2021' <- ssp

# Summer 2021
divmat_pa_year2 <- divmat_pa[which(dat_plot$Season.Year=="Summer 2021"),]
dat_plot_year2 <- dat_plot[which(dat_plot$Season.Year=="Summer 2021"),]
permmod <- adonis2(divmat_pa_year2 ~ Cover*Cropping.system, data = dat_plot_year2,  by = "term",
                   method="bray", permutations = 499, strata = dat_plot_year2$Replicate)
permmod2 <- round(as.data.frame(permmod), 3)
# append to a mod.effects dataframe
ssp <- paste0("F=", permmod2$F, ", Df=", permmod2$Df, ", p=", permmod2$`Pr(>F)`, ", R2=", permmod2$R2)
mod.effects$'Summer 2021' <- ssp

# Fall 2021
divmat_pa_year2 <- divmat_pa[which(dat_plot$Season.Year=="Fall 2021"),]
dat_plot_year2 <- dat_plot[which(dat_plot$Season.Year=="Fall 2021"),]
permmod <- adonis2(divmat_pa_year2 ~ Cover*Cropping.system, data = dat_plot_year2,  by = "term",
                   method="bray", permutations = 499, strata = dat_plot_year2$Replicate)
permmod2 <- round(as.data.frame(permmod), 3)
# append to a mod.effects dataframe
ssp <- paste0("F=", permmod2$F, ", Df=", permmod2$Df, ", p=", permmod2$`Pr(>F)`, ", R2=", permmod2$R2)
mod.effects$'Fall 2021' <- ssp

# Spring 2022
divmat_pa_year2 <- divmat_pa[which(dat_plot$Season.Year=="Spring 2022"),]
dat_plot_year2 <- dat_plot[which(dat_plot$Season.Year=="Spring 2022"),]
permmod <- adonis2(divmat_pa_year2 ~ Cover*Cropping.system, data = dat_plot_year2, by = "term",
                   method="bray", permutations = 499, strata = dat_plot_year2$Replicate)
permmod2 <- round(as.data.frame(permmod), 3)
# append to a mod.effects dataframe
ssp <- paste0("F=", permmod2$F, ", Df=", permmod2$Df, ", p=", permmod2$`Pr(>F)`, ", R2=", permmod2$R2)
mod.effects$'Spring 2022' <- ssp

# Summer 2022
divmat_pa_year2 <- divmat_pa[which(dat_plot$Season.Year=="Summer 2022"),]
dat_plot_year2 <- dat_plot[which(dat_plot$Season.Year=="Summer 2022"),]
permmod <- adonis2(divmat_pa_year2 ~ Cover*Cropping.system, data = dat_plot_year2,  by = "term",
                   method="bray", permutations = 499, strata = dat_plot_year2$Replicate)
permmod2 <- round(as.data.frame(permmod), 3)
# append to a mod.effects dataframe
ssp <- paste0("F=", permmod2$F, ", Df=", permmod2$Df, ", p=", permmod2$`Pr(>F)`, ", R2=", permmod2$R2)
mod.effects$'Summer 2022' <- ssp

# Fall 2022
divmat_pa_year2 <- divmat_pa[which(dat_plot$Season.Year=="Fall 2022"),]
dat_plot_year2 <- dat_plot[which(dat_plot$Season.Year=="Fall 2022"),]
permmod <- adonis2(divmat_pa_year2 ~ Cover*Cropping.system, data = dat_plot_year2,  by = "term",
                   method="bray", permutations = 499, strata = dat_plot_year2$Replicate)
permmod2 <- round(as.data.frame(permmod), 3)
# append to a mod.effects dataframe
ssp <- paste0("F=", permmod2$F, ", Df=", permmod2$Df, ", p=", permmod2$`Pr(>F)`, ", R2=", permmod2$R2)
mod.effects$'Fall 2022' <- ssp

# Spring 2023
divmat_pa_year2 <- divmat_pa[which(dat_plot$Season.Year=="Spring 2023"),]
dat_plot_year2 <- dat_plot[which(dat_plot$Season.Year=="Spring 2023"),]
permmod <- adonis2(divmat_pa_year2 ~ Cover*Cropping.system, data = dat_plot_year2, by = "term",
                   method="bray", permutations = 499, strata = dat_plot_year2$Replicate)
permmod2 <- round(as.data.frame(permmod), 3)
# append to a mod.effects dataframe
ssp <- paste0("F=", permmod2$F, ", Df=", permmod2$Df, ", p=", permmod2$`Pr(>F)`, ", R2=", permmod2$R2)
mod.effects$'Spring 2023' <- ssp

# Summer 2023
divmat_pa_year2 <- divmat_pa[which(dat_plot$Season.Year=="Summer 2023"),]
dat_plot_year2 <- dat_plot[which(dat_plot$Season.Year=="Summer 2023"),]
permmod <- adonis2(divmat_pa_year2 ~ Cover*Cropping.system, data = dat_plot_year2,  by = "term",
                   method="bray", permutations = 499, strata = dat_plot_year2$Replicate)
permmod2 <- round(as.data.frame(permmod), 3)
# append to a mod.effects dataframe
ssp <- paste0("F=", permmod2$F, ", Df=", permmod2$Df, ", p=", permmod2$`Pr(>F)`, ", R2=", permmod2$R2)
mod.effects$'Summer 2023' <- ssp

# Fall 2023
divmat_pa_year2 <- divmat_pa[which(dat_plot$Season.Year=="Fall 2023"),]
dat_plot_year2 <- dat_plot[which(dat_plot$Season.Year=="Fall 2023"),]
permmod <- adonis2(divmat_pa_year2 ~ Cover*Cropping.system, data = dat_plot_year2,  by = "term",
                   method="bray", permutations = 499, strata = dat_plot_year2$Replicate)
permmod2 <- round(as.data.frame(permmod), 3)
# append to a mod.effects dataframe
ssp <- paste0("F=", permmod2$F, ", Df=", permmod2$Df, ", p=", permmod2$`Pr(>F)`, ", R2=", permmod2$R2)
mod.effects$'Fall 2023' <- ssp

# Spring 2024
divmat_pa_year2 <- divmat_pa[which(dat_plot$Season.Year=="Spring 2024"),]
dat_plot_year2 <- dat_plot[which(dat_plot$Season.Year=="Spring 2024"),]
permmod <- adonis2(divmat_pa_year2 ~ Cover*Cropping.system, data = dat_plot_year2, by = "term",
                   method="bray", permutations = 499, strata = dat_plot_year2$Replicate)
permmod2 <- round(as.data.frame(permmod), 3)
# append to a mod.effects dataframe
ssp <- paste0("F=", permmod2$F, ", Df=", permmod2$Df, ", p=", permmod2$`Pr(>F)`, ", R2=", permmod2$R2)
mod.effects$'Spring 2024' <- ssp

# Summer 2024
divmat_pa_year2 <- divmat_pa[which(dat_plot$Season.Year=="Summer 2024"),]
dat_plot_year2 <- dat_plot[which(dat_plot$Season.Year=="Summer 2024"),]
permmod <- adonis2(divmat_pa_year2 ~ Cover*Cropping.system, data = dat_plot_year2,  by = "term",
                   method="bray", permutations = 499, strata = dat_plot_year2$Replicate)
permmod2 <- round(as.data.frame(permmod), 3)
# append to a mod.effects dataframe
ssp <- paste0("F=", permmod2$F, ", Df=", permmod2$Df, ", p=", permmod2$`Pr(>F)`, ", R2=", permmod2$R2)
mod.effects$'Summer 2024' <- ssp

# Fall 2024
divmat_pa_year2 <- divmat_pa[which(dat_plot$Season.Year=="Fall 2024"),]
dat_plot_year2 <- dat_plot[which(dat_plot$Season.Year=="Fall 2024"),]
permmod <- adonis2(divmat_pa_year2 ~ Cover*Cropping.system, data = dat_plot_year2,  by = "term",
                   method="bray", permutations = 499, strata = dat_plot_year2$Replicate)
permmod2 <- round(as.data.frame(permmod), 3)
# append to a mod.effects dataframe
ssp <- paste0("F=", permmod2$F, ", Df=", permmod2$Df, ", p=", permmod2$`Pr(>F)`, ", R2=", permmod2$R2)
mod.effects$'Fall 2024' <- ssp


# export model results
write.csv(mod.effects, paste0("Model-output/permanova/",name,"_permanova_by season.csv"))


mod.effects_sig <- mod.effects[1:3,]
for(i in 2:ncol(mod.effects_sig)){
  df <- strsplit(mod.effects_sig[,i], ", ", fixed=T)
  df2 <- c(df[[1]][3], df[[2]][3], df[[3]][3])
  rs <- c(df[[1]][4], df[[2]][4], df[[3]][4])
  df3 <- strsplit(df2, "p=", fixed=T)
  rs3 <- strsplit(rs, "R2=", fixed=T)
  df4 <- as.numeric(c(df3[[1]][2], df3[[2]][2], df3[[3]][2]))
  rs4 <- round(as.numeric(c(rs3[[1]][2], rs3[[2]][2], rs3[[3]][2]))*100,2)
  sig <- rep("NA", length(df4))
  for(t in 1:length(df4)){
    if(df4[t] < 0.05){sig[t] <- "*"}
    if(df4[t] < 0.01){sig[t] <- "**"}
    if(df4[t] < 0.001){sig[t] <- "***"}
  }
  mod.effects_sig[,i] <- paste(rs4, sig)
}

write.csv(mod.effects_sig, paste0("Model-output/permanova/",name,"_permanova_by season_*.csv"))








####################################
### Soil effects on communities: RDA
####################################

# import data for abundances and taxonomy
taxabund <- read.csv(paste0("Raw-data/sequence/",sampdate,"/",name,"/seqtab/6seqtab.csv"))
# import data for diversity
divdat <- read.csv(paste0("Processed-data/sequences/", region, "-diversity.csv"))
# import data for diversity
divmat_pa <- read.csv(paste0("Raw-data/sequence/",sampdate,"/divmat_pa_",region,".csv"), row.names = 1)


divdat$Season.Year <- as.factor(divdat$Season.Year)
divdat$Season.Year <- factor(divdat$Season.Year, levels(divdat$Season.Year)[c(1,6,10,2,7,11,3,8,12,4,9,13,5)])
divdat$Season <- as.factor(divdat$Season )
divdat$Season <- factor(divdat$Season, levels(divdat$Season)[c(2,3, 1)])
divdat$Cover <- as.factor(divdat$Cover)
divdat$Cover <- factor(divdat$Cover, levels(divdat$Cover)[c(2,1,4,5,3)])
divdat$Cropping.system <- as.factor(divdat$Cropping.system)
divdat$Cropping.system <- factor(divdat$Cropping.system, levels(divdat$Cropping.system)[c(1,4,3,2)])
divdat$Replicate <- as.factor(divdat$Replicate)
divdat$Year <- as.factor(divdat$Year)
divdat <- divdat[which(divdat$new.name %in% rownames(divmat_pa)),]


soildat$Sample <- paste0(soildat$Plot, "_", soildat$Season, "_", soildat$Year)
soildat2 <- soildat[which(soildat$Sample %in% rownames(divmat_pa)),]
rownames(soildat2) <- soildat2$Sample
soildat2 <- soildat2[rownames(divmat_pa),]


###### rda: visualize how communities vary with treatment and soil properties 

# first determine which properties are co-linear so that we 
# can exclude them from the rda if needed
res <- cor(soildat2[,c("GMC", "noN", "nhN", "inorgN", "WEC", "MBC", "BG",  "NAG", "PHOS")], 
           method="pearson", use="pairwise.complete.obs")
res2 <- as.data.frame(round(res, 2))
#install.packages("corrplot") 
library(corrplot)
pdf("Figures/soil correlations/soil-property-correlations.pdf", width=7, height=7)
corrplot(res, type = "upper", order = "hclust", method = "number", diag = FALSE, 
         tl.col = "black", tl.srt = 45)
dev.off()
# inspect plot



### do soil properties and treatments influence microbial community composition?
mod.effects <- data.frame(effect=c("model", "residual")) 



# overall rda
keep <- which(complete.cases(soildat2[,c("GMC", "noN", "nhN", "WEC", "MBC", "BG",  "NAG", "PHOS")]))
# myrdaall <- capscale(formula = divmat_pa[keep,] ~ GMC + noN + nhN + WEC + MBC + BG + NAG + PHOS,
#                      distance = "bray",
#                      data = soildat2[keep, ], na.action = "na.omit")
# myrdaall_s <- ordistep(myrdaall, perm.max = 10)

formula_s <- str_split((summary(myrdaall_s))[5], pattern="~")
formula_s <- str_split(formula_s[[1]][2], pattern=",")
preds <- formula_s[[1]][1]
myrdaall <- capscale(formula = as.formula(paste("divmat_pa[keep,] ~ ", preds)),
                     distance = "bray", 
                     data = soildat2[keep, ], na.action = "na.omit")
vif.cca(myrdaall)
mod <- round(as.data.frame(anova(myrdaall)), 3)
ssp <- paste0("F=", mod$F, ", Df=", mod$Df, ", p=", mod$`Pr(>F)`, ", R2=", round(mod$SumOfSqs/sum(mod$SumOfSqs[1:2]),3))

# append to a mod.effects dataframe
mod.effects$myrdaall <- ssp

## rda all
# factors associated with each axis
write.csv(as.data.frame(summary(myrdaall)$biplot), paste0("Model-output/db-rda/",name,"_all_loadings.csv"))
# site coordinates
write.csv(as.data.frame(summary(myrdaall)$site), paste0("Model-output/db-rda/",name,"_all_axes.csv"))
# variation explained by first two axes
write.csv(summary(myrdaall)$cont$importance, paste0("Model-output/db-rda/",name,"_all_variance-explained.csv"))
# variation explained by each soil property
ef <- envfit(myrdaall, na.omit(soildat2[, c("GMC","noN", "nhN", "WEC", "MBC", "BG", "NAG", "PHOS")]), choices = c(1,2), permutations = 0, na.rm=TRUE)
write.csv(round(cbind(ef$vectors$arrows, ef$vectors$r),2), paste0("Model-output/db-rda/",name,"_all_variance-by-predictors.csv"))




#### visualize: all by cropping system
# with soil property eigenvectors
rdadat <- as.data.frame(summary(myrdaall)$sites) # sites
rem <- which(complete.cases(soildat2[, c("GMC", "noN", "nhN", "WEC", "MBC", "BG",  "NAG", "PHOS")])==FALSE)
rdadat <- cbind(rdadat, divdat[-rem, c("Season.Year", "Cover", "Cropping.system", "Year", "Season")]) #
df2  <- data.frame(summary(myrdaall)$biplot) # loadings: only keep top rated for first two axes
df3 <- df2
var1 <- round(read.csv(paste0("Model-output/db-rda/",name,"_all_variance-explained.csv"), row.names = 1)[2,1]*100, 2)
var2 <- round(read.csv(paste0("Model-output/db-rda/",name,"_all_variance-explained.csv"), row.names = 1)[2,2]*100, 2)

# save.image("Model-output/db-rda/Bacteria-all-dat.RData")
load("Model-output/db-rda/Bacteria-all-dat.RData") # start from here if revising plot

# Extract species (bacterial families) scores
species_scores <- scores(myrdaall, display = "species")

# Convert to data frame for ggplot
df_bact <- as.data.frame(species_scores[, 1:2])  # Only take first 2 axes (RDA1, RDA2)

# Optionally scale for visual clarity (e.g., shrink arrows to reduce clutter)
scaling_factor <- 200
df_bact$CAP1 <- df_bact[, 1] * scaling_factor
df_bact$CAP2 <- df_bact[, 2] * scaling_factor
df_bact$label <- taxabund$Phylum # could try with other biological levels
df_bact <- df_bact %>%
  group_by(label) %>%
  summarise(CAP1 = mean(CAP1, na.rm = TRUE),
            CAP2 = mean(CAP2, na.rm = TRUE)) %>% 
  mutate(sum = abs(CAP1)+(CAP2)) %>% 
  arrange(desc(sum)) %>% 
  top_n(10)


pdf(paste0("Figures/microbial-diversity/rda/",name,"_dbrda-all1.pdf"),  height=5.5, width=7)
ggplot(rdadat, aes(CAP1, CAP2)) +
  geom_point(aes(fill=Season.Year, shape=Cropping.system)) +
  theme_bw() +
  scale_fill_manual(values = ssnyrcols) +
  scale_shape_manual(values = c(21:24)) +
  guides(fill=guide_legend(title="Season and year", override.aes=list(shape=21)),
         shape=guide_legend(title="Cropping system")) +
  labs(x=paste0("RDA1 (", var1, "%)"), y=paste0("RDA2 (", var2, "%)"), title=name) + # amend according to variance explained
  lims(x = c(-2,2), y = c(-2,3)) +
  # Bacterial family vectors
  geom_segment(data = df_bact, aes(x = 0, xend = CAP1, y = 0, yend = CAP2), 
               color = "magenta3", size = 0.75,
               arrow = arrow(length = unit(0.02, "npc"))) +
  geom_label(data = df_bact, aes(x = CAP1, y = CAP2, label = label), 
             hjust = 0.5 * (1 - sign(df_bact$CAP1)),
             vjust = 0.5 * (1 - sign(df_bact$CAP2)),
             color = "magenta3", label.size = NA, fill = alpha(c("white"),0.25), size = 2) +
  geom_segment(data=df3, aes(x=0, xend=CAP1, y=0, yend=CAP2), color="white", size=1.25,arrow=arrow(length=unit(0.02,"npc"))) +
  geom_segment(data=df3, aes(x=0, xend=CAP1, y=0, yend=CAP2), color="black", size=1,arrow=arrow(length=unit(0.02,"npc"))) +
  geom_label(data=df3, aes(x=df3$CAP1,y=df3$CAP2,label=rownames(df3)), 
             hjust=0.5*(1-sign(df3$CAP1)),vjust=0.5*(1-sign(df3$CAP2)), fill = "white", size=2) 
  

dev.off()


#### visualize: all by cover
# with soil property eigenvectors
pdf(paste0("Figures/microbial-diversity/rda/",name,"_dbrda-all3.pdf"), height=5.5, width=7)
ggplot(rdadat, aes(CAP1, CAP2)) +
  geom_point(aes(fill=Season.Year, shape=Cover)) +
  theme_bw() +
  scale_fill_manual(values = ssnyrcols) +
  scale_shape_manual(values = c(21:25)) +
  guides(fill=guide_legend(title="Season and year", override.aes=list(shape=21)),
         shape=guide_legend(title="Cover")) +
  labs(x=paste0("RDA1 (", var1, "%)"), y=paste0("RDA2 (", var2, "%)"), title=name) + # amend according to variance explained
  lims(x = c(-2,2), y = c(-2,3)) +
  # Bacterial family vectors
  geom_segment(data = df_bact, aes(x = 0, xend = CAP1, y = 0, yend = CAP2), 
               color = "magenta3", size = 0.75,
               arrow = arrow(length = unit(0.02, "npc"))) +
  geom_label(data = df_bact, aes(x = CAP1, y = CAP2, label = label), 
             hjust = 0.5 * (1 - sign(df_bact$CAP1)),
             vjust = 0.5 * (1 - sign(df_bact$CAP2)),
             color = "magenta3", label.size = NA, fill = alpha(c("white"),0.25), size = 2) +
  geom_segment(data=df3, aes(x=0, xend=CAP1, y=0, yend=CAP2), color="white", size=1.25,arrow=arrow(length=unit(0.02,"npc"))) +
  geom_segment(data=df3, aes(x=0, xend=CAP1, y=0, yend=CAP2), color="black", size=1,arrow=arrow(length=unit(0.02,"npc"))) +
  geom_label(data=df3, aes(x=df3$CAP1,y=df3$CAP2,label=rownames(df3)), 
             hjust=0.5*(1-sign(df3$CAP1)),vjust=0.5*(1-sign(df3$CAP2)), fill = "white", size=2) 
dev.off()


# export mod effects
write.csv(t(mod.effects), paste0("Model-output/db-rda/",name,"_dbrda_all timepoints.csv"))












####################################
### Soil property influence on communities for each seasonXyear
####################################

# import data for diversity
divdat <- read.csv(paste0("Processed-data/sequences/", region,"-diversity.csv"))
divdat$Season.Year <- as.factor(divdat$Season.Year)
divdat$Season <- as.factor(divdat$Season )
divdat$Season <- factor(divdat$Season, levels(divdat$Season)[c(2,3, 1)])
divdat$Cover <- as.factor(divdat$Cover)
divdat$Cover <- factor(divdat$Cover, levels(divdat$Cover)[c(2,1,4,5,3)])
divdat$Cropping.system <- as.factor(divdat$Cropping.system)
divdat$Cropping.system <- factor(divdat$Cropping.system, levels(divdat$Cropping.system)[c(1,4,3,2)])
divdat$Replicate <- as.factor(divdat$Replicate)
divdat$Year <- as.factor(divdat$Year)
divdat <- divdat[which(divdat$new.name %in% rownames(divmat_pa)),]
# import data for diversity
divmat_pa0 <- read.csv(paste0("Raw-data/sequence/",sampdate,"/divmat_pa_", region, ".csv"), row.names = 1)
# soil data
soildat$Sample <- paste0(soildat$Plot, "_", soildat$Season, "_", soildat$Year)

resp_list <- list()




### RDA: do soil properties and treatments influence microbial community composition?
# create data frame for results
mod.effects <- data.frame(effect=c("model", "residual")) 
# timepoints
seasons <- levels(soildat$season.year)
# variance inflation factors
vifs <- list(NA)
# significance of overall model
resp_list <- list(NA)
times <- levels(soildat$season.year)
time_l <- length(times)
for(t in 1:time_l){
  new_list <- list(NA)
  resp_list <- c(resp_list, new_list)
}
resp_list <- resp_list[-1]
names(resp_list) <- times
# create lists to plot rdas later
df3 <-  vector('list', length(seasons))
rdadat <-  vector('list', length(seasons))
myrdaall <- vector('list', length(seasons))
divdat0 <- vector('list', length(seasons))
var1 <- rep(NA, length(seasons))
var2 <- rep(NA, length(seasons))

for(i in 1:length(seasons)) {  
  
  time <- seasons[i]
  ti_me <- gsub(pattern = " ", replacement = "_", time)
  
  divdat0[[i]] <- divdat[which(divdat$Season.Year == time),]
  divmat_pa <- divmat_pa0[which(str_detect(rownames(divmat_pa0), pattern=ti_me)),]
  
  soildat2 <- soildat[which(soildat$Sample %in% rownames(divmat_pa)),]
  rownames(soildat2) <- soildat2$Sample
  soildat2 <- soildat2[rownames(divmat_pa),]
  
  print(paste0("Starting ", time))
  
  ###### rda: visualize how communities vary with treatment and soil properties
  
  # first determine which properties are co-linear so that we 
  # can exclude them from the rda if needed
  res <- cor(soildat2[which(soildat2$season.year==time),c("GMC", "noN", "nhN", "inorgN", "WEC", "MBC", "BG",  "NAG", "PHOS")], 
             method="pearson", use="pairwise.complete.obs")
  res2 <- as.data.frame(round(res, 2))
  #install.packages("corrplot") 
  library(corrplot)
  pdf(paste0("Figures/soil correlations/soil-property-correlations",time,".pdf"), width=7, height=7)
  corrplot(res, type = "upper", order = "hclust", method = "number", diag = FALSE, 
           tl.col = "black", tl.srt = 45)
  dev.off()
  # inspect plot
  
  
  # overall rda
  keep <- which(complete.cases(soildat2[,c("GMC", "noN", "nhN", "WEC", "MBC", "BG",  "NAG", "PHOS")]))
  myrdaall[[i]] <- capscale(formula = divmat_pa[keep,] ~ GMC + noN + nhN + WEC + MBC + BG + NAG + PHOS,
                       distance = "bray", 
                       data = soildat2[keep, ], na.action = "na.omit")
  myrdaall_s <- ordistep(myrdaall[[i]], perm.max = 10) 
  formula_s <- str_split((summary(myrdaall_s))[5], pattern="~")
  formula_s <- str_split(formula_s[[1]][2], pattern=",")
  preds <- formula_s[[1]][1]
  myrdaall[[i]] <- capscale(formula = as.formula(paste("divmat_pa[keep,] ~ ", preds)),
                       distance = "bray", 
                       data = soildat2[keep, ], na.action = "na.omit")
  vifs_time <- list(vif.cca(myrdaall[[i]]))
  vifs <- c(vifs, vifs_time)
  mod <- round(as.data.frame(anova(myrdaall[[i]])), 3)
  ps <- mod$`Pr(>F)`[1]
  if(ps >= 0.05){resp_list[[i]] <- "N.S."}
  else if(ps >= 0.01){resp_list[[i]] <- "*"}
  else if(ps >= 0.001){resp_list[[i]] <- "**"}
  else {resp_list[[i]] <- "***"}
  ssp <- paste0("F=", mod$F, ", Df=", mod$Df, ", p=", mod$`Pr(>F)`, ", R2=", round(mod$SumOfSqs/sum(mod$SumOfSqs[1:2]),3))
  
  # append to a mod.effects dataframe
  mod.effects[,paste0("myrda",ti_me)] <- ssp
  
  ## rda all
  # factors associated with each axis
  write.csv(as.data.frame(summary(myrdaall[[i]])$biplot), paste0("Model-output/db-rda/",name,"_all_loadings_", time, ".csv"))
  # site coordinates
  write.csv(as.data.frame(summary(myrdaall[[i]])$site), paste0("Model-output/db-rda/",name,"_all_axes_", time, ".csv"))
  # variation explained by first two axes
  write.csv(summary(myrdaall[[i]])$cont$importance, paste0("Model-output/db-rda/",name,"_all_variance-explained_", time, ".csv"))
  # variation explained by each soil property
  ef <- envfit(myrdaall[[i]], na.omit(soildat2[, c("GMC","noN", "nhN", "WEC", "MBC", "BG", "NAG", "PHOS")]), choices = c(1,2), permutations = 100, na.rm=TRUE)
  ef_df <- round(cbind(ef$vectors$arrows, ef$vectors$r, ef$vectors$pvals),4)
  colnames(ef_df)[3:4] <- c("r2", "pvalue")
  write.csv(ef_df, paste0("Model-output/db-rda/",name,"_all_variance-by-predictors_", time, ".csv"))
  
  
  #### visualize: all by cropping system and cover
  # with soil property eigenvectors
  rdadat0 <- as.data.frame(summary(myrdaall[[i]])$sites) # sites
  rownames(divdat0[[i]]) <- divdat0[[i]]$new.name
  rdadat[[i]] <- merge(divdat0[[i]], rdadat0, by = "row.names", all = TRUE)
  
  df2  <- data.frame(summary(myrdaall[[i]])$biplot) # loadings: only keep top rated for first two axes
  df3[[i]] <- df2
}

# save.image("Model-output/db-rda/Bacteria-dat.RData")
load("Model-output/db-rda/Bacteria-dat.RData") # start from here if revising figures

for(i in 1:length(df3)){
  time <- seasons[i]
  write.csv(df3[[i]], paste0("Model-output/db-rda/",name,"_all_model-predictors_", time, ".csv"))
}


# get variance explained
for(i in 1:length(seasons)){
  time <- seasons[i]
  var1[i] <- round(read.csv(paste0("Model-output/db-rda/",name,"_all_variance-explained_", time, ".csv"), row.names = 1)[3,1]*100, 2)
  var2[i] <- round(read.csv(paste0("Model-output/db-rda/",name,"_all_variance-explained_", time, ".csv"), row.names = 1)[3,2]*100, 2)
}


# rda plots
list_plots <- vector('list', length(seasons))
df_bact <- list()


# export multipanel rda figure
plotList <- lapply(
  1:length(seasons),
  function(key) {
    rdadat[[key]]$Cover <- as.factor(rdadat[[key]]$Cover)
    #rdadat[[key]]$Cover <- factor(rdadat[[key]]$Cover, levels(rdadat[[key]]$Cover)[c(2,4,1,5,3)])
    rdadat[[key]]$Cropping.system <- as.factor(rdadat[[key]]$Cropping.system)
    #rdadat[[key]]$Cropping.system <- factor(rdadat[[key]]$Cropping.system, levels(rdadat[[key]]$Cropping.system)[c(1,4,3,2)])
    # Extract species (bacterial families) scores
    species_scores <- scores(myrdaall[[key]], display = "species")
    
    # Convert to data frame for ggplot
    df_bact0 <- as.data.frame(species_scores[, 1:2])  # Only take first 2 axes (RDA1, RDA2)
    
    # Optionally scale for visual clarity (e.g., shrink arrows to reduce clutter)
    scaling_factor <- 200
    df_bact0$CAP1 <- df_bact0[, 1] * scaling_factor
    df_bact0$CAP2 <- df_bact0[, 2] * scaling_factor
    df_bact0$label <- taxabund$Phylum
    df_bact[[key]] <- df_bact0 %>%
      group_by(label) %>%
      summarise(CAP1 = mean(CAP1, na.rm = TRUE),
                CAP2 = mean(CAP2, na.rm = TRUE)) %>% 
      mutate(sum = abs(CAP1)+(CAP2)) %>% 
      arrange(desc(sum)) %>% 
      top_n(3)
    
    df3[[key]] <- df3[[key]] %>%
      mutate(
        quadrant = case_when(
          CAP1 > 0 & CAP2 > 0 ~ 45,
          CAP1 < 0 & CAP2 > 0 ~ -45,
          CAP1 < 0 & CAP2 < 0 ~ 45,
          CAP1 > 0 & CAP2 < 0 ~ -45
        )
      )
    
    # Need to assign the plot to a variable because 
    # you want to generate the plot AND save to file 
    x <- ggplot(rdadat[[key]], aes(CAP1, CAP2)) +
      geom_point(aes(fill=Cropping.system, shape=Cover), size=2, alpha = 0.5) +
      theme_bw() +
      scale_fill_manual(values = cropsyscols) +
      scale_shape_manual(values = c(21:25)) +
      guides(fill=guide_legend(title="Cropping system", override.aes=list(shape=21, size=3)),
             shape=guide_legend(title="Cover", override.aes = list(size=3))) +
      labs(x=paste0("RDA1 (", var1[key], "%)"), y=paste0("RDA2 (", var2[key], "%)"), title=seasons[key]) + # amend according to variance explained
      theme(legend.text = element_text(size = 14), 
            legend.title = element_text(size = 16)) +
      lims(x = c(-2,2), y = c(-2,3)) +
      # Fungal family vectors
      # geom_segment(data = df_bact[[key]], aes(x = 0, xend = CAP1, y = 0, yend = CAP2), 
      #              color = "magenta3", size = 0.75,
      #              arrow = arrow(length = unit(0.02, "npc"))) +
      # geom_label(data = df_bact[[key]], aes(x = CAP1, y = CAP2, label = label), 
      #            hjust = 0.5 * (1 - sign(df_bact[[key]]$CAP1)),
      #            vjust = 0.5 * (1 - sign(df_bact[[key]]$CAP2)),
      #            color = "magenta3", label.size = NA, fill = alpha(c("white"),0.25), size = 3, fontface="italic") +
      # Soil property vectors
      geom_segment(data=df3[[key]], aes(x=0, xend=CAP1, y=0, yend=CAP2), color="white", size=1,arrow=arrow(length=unit(0.02,"npc"))) +
      geom_segment(data=df3[[key]], aes(x=0, xend=CAP1, y=0, yend=CAP2), color="blue", size=0.75,arrow=arrow(length=unit(0.02,"npc"))) +
      geom_text(data=df3[[key]], aes(x=df3[[key]]$CAP1,# + 0.4*sign(df3[[key]]$CAP1)*df3[[key]]$CAP1, 
                                     y=df3[[key]]$CAP2,# + 0.4*sign(df3[[key]]$CAP2)*df3[[key]]$CAP2, 
                                     angle = df3[[key]]$quadrant), 
                label=rownames(df3[[key]]), size=3.2, fontface="bold", color = "white", 
                hjust=0.5*(1-sign(df3[[key]]$CAP1)),
                vjust=0.5*(1-sign(df3[[key]]$CAP2))) +
      geom_text(data=df3[[key]], aes(x=df3[[key]]$CAP1,# + 0.4*sign(df3[[key]]$CAP1)*df3[[key]]$CAP1, 
                                     y=df3[[key]]$CAP2,# + 0.4*sign(df3[[key]]$CAP2)*df3[[key]]$CAP2, 
                                     angle = df3[[key]]$quadrant), 
                label=rownames(df3[[key]]), size=3, fontface="bold", color = "blue", 
                hjust=0.5*(1-sign(df3[[key]]$CAP1)),
                vjust=0.5*(1-sign(df3[[key]]$CAP2)))
  }
)



allplots <- ggarrange(plotlist=plotList, 
                      labels = "AUTO", common.legend = T, legend = "right",
                      ncol = 3, nrow=5)

png(paste0("Figures/microbial-diversity/rda/",name,"_dbrda-all1_across-timepoints.png"),  height=4000, width=3500, res = 300)
allplots 
dev.off()


# examine vifs
vifs
vifs0 <- vifs[-1]
vifs2 <-  bind_rows(vifs0, .id = "timepoint")
vifs3 <- vifs2[,-1]
vifs3 <- as.data.frame(vifs3)
rownames(vifs3) <- seasons
vifs3

# examine total variance explained
names(resp_list) <- seasons
resp <- bind_cols(resp_list, .id = "timepoint")



# export results
write.csv(t(mod.effects), paste0("Model-output/db-rda/*",name,"_dbrda_by timepoint.csv"))
write.csv(vifs3, paste0("Model-output/db-rda/*",name,"_dbrda_by timepoint_vifs.csv"))
write.csv(resp, paste0("Model-output/db-rda/*",name,"_dbrda_by timepoint_total-r2-sig.csv"))





