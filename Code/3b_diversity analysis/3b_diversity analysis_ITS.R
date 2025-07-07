##################################################
###### 1.  diversity analysis




library(vegan)
#install.packages("devtools")
#devtools::install_github("jfq3/QsRutils", build_vignettes = TRUE)
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
library(MuMIn)


# select colors
cropsyscols <- c("darkgoldenrod4", "darkolivegreen", "darkolivegreen3", "deepskyblue3")
covercols <- c("chocolate", "goldenrod", "forestgreen", "darkcyan", "darkorchid4")
ssnyrcols <- c("white", 
               RColorBrewer::brewer.pal(n = 4, name = "Purples")[2:4], 
               RColorBrewer::brewer.pal(n = 4, name = "YlGn")[2:4], 
               RColorBrewer::brewer.pal(n = 4, name = "Reds")[2:4], 
               RColorBrewer::brewer.pal(n = 4, name = "Blues")[2:4])
yrcols <- c("white", "purple", "aquamarine3", "brown2", "blue")



nplots <- 80 # number of plots within the field site




## there are 13 sampling points:
# fall 2020, (baseline)
# spring 2021, summer 2021, fall 2021, 
# spring 2022, summer 2022, fall 2022
# spring 2023, summer 2023, fall 2023
# spring 2024, summer 2024, fall 2024


################  Soil data  ################  

# field plot metadata
plotdat <- read.csv("Raw-data/WTREC-CDP_plot-data.csv")
# field plot soil property data
soildat <- read.csv("Processed-data/WTREC-CDP_processed-data-simple.csv")
# subset soildat to 0-10cm depth
soildat <- soildat[which(soildat$Depth=="0-10 cm"),] 
soildat$season.year <- paste(soildat$Season, soildat$Year, sep=" ")
soildat$season.year <- as.factor(soildat$season.year )
soildat$season.year <- factor(soildat$season.year, levels(soildat$season.year)[c(1,6,10,2,7,11,3,8,12,4,9,13,5)])






################  Fungal diversity ################  

sampdate <- c("All_agmicrobiome") 
region <- "ITS"
name <- "Fungi"
nm <- "fun"
# number of sampling periods for which we have sequence data
nsamps <- 13


divdat <- read.csv(paste0("Processed-data/sequences/", region, "-diversity.csv")) 



# merge datasets 
soildat$new.name <- paste(soildat$Plot, soildat$Season, soildat$Year, sep="_") 
unique(divdat$new.name == soildat$new.name) # make sure soil and diversity data are ordered in the same way

## add diversity data to soil data
soildat[,which(substr(colnames(soildat), 1, 3) == nm)] <- divdat[, c("richness", "shannon.div", "simpson.div", "invsimpson.div", "evenness")]
soildat$amf_relabund <- divdat$amf_relabund
soildat$plantpath_relabund <- divdat$plantpath_relabund
soildat$sap_relabund <- divdat$sap_relabund

# export merged data 
write.csv(soildat, paste0("Processed-data/sequences/soil_", region, "-diversity.csv"))



##################### MODELING APPROACH #####################
### We have a split-plot block design, 
### with 4 main plot levels (cropping system)
### and 5 subplot levels (cover)  
### with each treatment combo represented within each of 4 blocks (Replicate)
### and sampled at multiple timepoints (Season.Year)
### There are three predictor variables:
### 1. Replicate (4 levels), random effect
### 2. Cropping system (4 levels), main plot
### 3. Cover (5 levels), subplot
### and we will model separately for each Year or Season.Year


responses <- c("fun_richness", "fun_shannon.div", "fun_simpson.div", "fun_invsimpson.div", "fun_evenness", "amf_relabund", "plantpath_relabund", "sap_relabund")
resp_l <- length(responses)

times <- levels(soildat$season.year)
time_l <- length(times)


dat2 <- soildat


################################################

# dataframe for model output
resp <- data.frame(resp = rep(responses),
                   normality_pvalue = rep(NA, resp_l),
                   normality_W = rep(NA, resp_l),
                   cover = rep(NA, resp_l),
                   cropsys = rep(NA, resp_l),
                   time = rep(NA, resp_l),
                   cover.cropsys = rep(NA, resp_l),
                   cover.time = rep(NA, resp_l),
                   cropsys.time = rep(NA, resp_l),
                   cover.cropsys.time = rep(NA, resp_l),
                   rsq = rep(NA, resp_l),
                  transformed = rep(NA, resp_l),
                  cover.p = rep(NA, resp_l),
                  cropsys.p = rep(NA, resp_l),
                  time.p = rep(NA, resp_l),
                  cover.cropsys.p = rep(NA, resp_l),
                  cover.time.p = rep(NA, resp_l),
                  cropsys.time.p = rep(NA, resp_l),
                  cover.cropsys.time.p = rep(NA, resp_l))
resp_list <- resp


# overall model
for (i in 1:resp_l) { 
  print(paste("Starting ", resp$resp[i], " (", i, " of ", resp_l,")"))
  response_var <- resp$resp[i]  # Current response variable
  
  fit <- lme(
    fixed = as.formula(paste("(", response_var, ") ~ Cover * Cropping.system * season.year")),
    random = ~ 1 | Rep,
    #weights = varIdent(form = ~ 1 | season.year),
    data = dat2,
    na.action = na.omit)
  
  # check for normality
  resid <- residuals(fit)
  resp_list$normality_pvalue[i] <- shapiro.test(resid)$p.value
  resp_list$normality_W[i] <- shapiro.test(resid)$statistic
  resp_list$transformed[i] <- "no"
  
  if(resp_list$normality_pvalue[i] < 0.05){
    # Fit the model
    fit <- lme(
      fixed = as.formula(paste("log(", response_var, "+1) ~ Cover * Cropping.system * season.year")),
      random = ~ 1 | Rep,
      #weights = varIdent(form = ~ 1 | season.year),
      data = dat2,
      na.action = na.omit)
    # check for normality
    resid <- residuals(fit)
    resp_list$normality_pvalue[i] <- shapiro.test(resid)$p.value
    resp_list$normality_W[i] <- shapiro.test(resid)$statistic
    resp_list$transformed[i] <- "yes"
  }
  
  # export model stats
  mod.out <- as.data.frame(anova(fit, type="marginal")) 
  mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
  resp_list[i,c(4:10)] <- paste0(round(mod.out2$`F-value`,2))[2:8]
  resp_list$rsq[i] <- paste0(round(r.squaredGLMM(fit)[1],2), ", ", round(r.squaredGLMM(fit)[2],2)) # marginal, conditional
  
  # replace p's with asterisks to denote significance
  ps <- mod.out2$`p-value`[2:8]
  for (k in 1:7){
    if(ps[k] >= 0.05){resp_list[i, k+12] <- ""}
    else if(ps[k] >= 0.01){resp_list[i, k+12] <- "*"}
    else if(ps[k] >= 0.001){resp_list[i, k+12] <- "**"}
    else {resp_list[i, k+12] <- "***"}
  }
  
  ### tukey: Cropping.system * Timepoint
  test1 <- emmeans(fit, ~ Cropping.system*season.year)
  testlet <- cld(test1, type = "response", Letters = "ABCDEFGH", reversed = TRUE)
  testlet <- testlet %>% dplyr::mutate_each_(funs(factor(.)),c("Cropping.system","season.year"))
  testlet$Cropping.system <- factor(testlet$Cropping.system, levels(testlet$Cropping.system)[c(1,4,3,2)])
  #testlet$Cover <- factor(testlet$Cover, levels(testlet$Cover)[c(2,4,1,5,3)])
  testlet <- testlet[order(testlet$season.year, testlet$Cropping.system),]
  testlet$.group <- gsub(" ", "", testlet$.group)
  if(length(unique(testlet$.group))==1) {testlet$.group <- rep(" ", length(testlet$.group))}
  write.csv(testlet, paste0("Model-output/anova_", region, "/cropsys by season.year/",response_var, ".csv"))
  
  ### tukey: Cover * Timepoint
  test1 <- emmeans(fit, ~ Cover*season.year)
  testlet <- cld(test1, type = "response", Letters = "ABCDEFGH", reversed = TRUE)
  testlet <- testlet %>% dplyr::mutate_each_(funs(factor(.)),c("Cover","season.year"))
  #testlet$season.year <- factor(testlet$season.year, levels(testlet$season.year)[c(1,4,3,2)])
  testlet$Cover <- factor(testlet$Cover, levels(testlet$Cover)[c(2,4,1,5,3)])
  testlet <- testlet[order(testlet$season.year, testlet$Cover),]
  testlet$.group <- gsub(" ", "", testlet$.group)
  if(length(unique(testlet$.group))==1) {testlet$.group <- rep(" ", length(testlet$.group))}
  write.csv(testlet, paste0("Model-output/anova_", region, "/cover by season.year/",response_var, ".csv"))
  
  ### tukey: Cover * Cropping system
  test1 <- emmeans(fit, ~ Cover*Cropping.system)
  testlet <- cld(test1, type = "response", Letters = "ABCDEFGH", reversed = TRUE)
  testlet <- testlet %>% dplyr::mutate_each_(funs(factor(.)),c("Cover","Cropping.system"))
  testlet$Cropping.system <- factor(testlet$Cropping.system, levels(testlet$Cropping.system)[c(1,4,3,2)])
  testlet$Cover <- factor(testlet$Cover, levels(testlet$Cover)[c(2,4,1,5,3)])
  testlet <- testlet[order(testlet$Cropping.system, testlet$Cover),]
  testlet$.group <- gsub(" ", "", testlet$.group)
  if(length(unique(testlet$.group))==1) {testlet$.group <- rep(" ", length(testlet$.group))}
  write.csv(testlet, paste0("Model-output/anova_", region, "/cover by cropsys/",response_var, ".csv"))
  
  ### tukey: Timepoint 
  test1 <- emmeans(fit, ~ season.year)
  testlet <- cld(test1, type = "response", Letters = "ABCDEFG", reversed = TRUE)
  testlet <- testlet %>% dplyr::mutate_each_(funs(factor(.)),c("season.year"))
  #testlet$season.year <- factor(testlet$season.year, levels(testlet$season.year)[c(2,4,1,5,3)])
  testlet <- testlet[order(testlet$season.year),]
  testlet$.group <- gsub(" ", "", testlet$.group)
  if(length(unique(testlet$.group))==1) {testlet$.group <- rep(" ", length(testlet$.group))}
  write.csv(testlet, paste0("Model-output/anova_", region, "/season.year/",response_var, ".csv"))
  
  ### tukey: Cropping.system 
  test1 <- emmeans(fit, ~ Cropping.system)
  testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
  testlet <- testlet %>% dplyr::mutate_each_(funs(factor(.)),c("Cropping.system"))
  testlet$Cropping.system <- factor(testlet$Cropping.system, levels(testlet$Cropping.system)[c(1,4,3,2)])
  testlet <- testlet[order(testlet$Cropping.system),]
  testlet$.group <- gsub(" ", "", testlet$.group)
  if(length(unique(testlet$.group))==1) {testlet$.group <- rep(" ", length(testlet$.group))}
  write.csv(testlet, paste0("Model-output/anova_", region, "/cropsys/",response_var,".csv"))
  
  ### tukey: Cover 
  test1 <- emmeans(fit, ~ Cover)
  testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
  testlet <- testlet %>% dplyr::mutate_each_(funs(factor(.)),c("Cover"))
  testlet$Cover <- factor(testlet$Cover, levels(testlet$Cover)[c(2,4,1,5,3)])
  testlet <- testlet[order(testlet$Cover),]
  testlet$.group <- gsub(" ", "", testlet$.group)
  if(length(unique(testlet$.group))==1) {testlet$.group <- rep(" ", length(testlet$.group))}
  write.csv(testlet, paste0("Model-output/anova_", region, "/cover/",response_var, ".csv"))
}

View(resp_list) # check
resp <- resp_list[order(resp_list$resp),]
write.csv(resp, paste0("Model-output/anova_", region, "/all-anova-overall.csv"))








################## model for each timepoint
responses <- c("fun_richness", "fun_shannon.div", "fun_simpson.div", "fun_invsimpson.div", "fun_evenness", "amf_relabund", "plantpath_relabund", "sap_relabund")
resp_l <- length(responses)

times <- levels(soildat$season.year)
time_l <- length(times)


dat2 <- soildat

# # Our model for each response variable will look like this:
# # with time as the repeated measure
# fit <- lme(resp ~ Cover * Cropping.system, # Fixed effects
#            random = ~ 1 | Rep, # Random intercept for Block (Rep)
#            weights = varIdent(form = ~ 1 | Cropping.system), # Allow unequal variances across treatments
#            data = dat2_time, # Dataset
#            na.action = na.omit # Handle missing values
# )



################################################

# dataframe for model output
resp <- data.frame(resp = rep(responses),
                   normality_pvalue = rep(NA, resp_l),
                   normality_W = rep(NA, resp_l),
                   cover = rep(NA, resp_l),
                   cropsys = rep(NA, resp_l),
                   cover.cropsys = rep(NA, resp_l),
                   rsq = rep(NA, resp_l),
                   transformed = rep(NA, resp_l),
                   cover_sig = rep(NA, resp_l),
                   cropsys_sig = rep(NA, resp_l),
                   cover.cropsys_sig = rep(NA, resp_l))
resp_list <- list(NA)
for(t in 1:time_l){
  new_list <- list(name=resp)
  resp_list <- c(resp_list, new_list)
}
resp_list <- resp_list[-1]
names(resp_list) <- times
# Loop through responses and fit models
for (i in 1:resp_l) { 
  response_var <- resp$resp[i]  # Current response variable
  
  for(t in 1:time_l){
    timepoint <- names(resp_list)[[t]]
    dat2_time <- dat2[which(dat2$season.year==timepoint),]  # Current timepoint
    
    # Fit the model
    fit <- lme(
      fixed = as.formula(paste("(", response_var, ") ~ Cover * Cropping.system")),
      random = ~ 1 | Rep,
      weights = varIdent(form = ~ 1 | Cropping.system),
      data = dat2_time,
      na.action = na.omit)
    # check for normality
    resid <- residuals(fit)
    resp_list[[t]]$normality_pvalue[i] <- shapiro.test(resid)$p.value
    resp_list[[t]]$normality_W[i] <- shapiro.test(resid)$statistic
    resp_list[[t]]$transformed[i] <- "no"
    
    if(resp_list[[t]]$normality_pvalue[i] < 0.05){
      # Fit the model
      fit <- lme(
        fixed = as.formula(paste("log(", response_var, "+1) ~ Cover * Cropping.system")),
        random = ~ 1 | Rep,
        weights = varIdent(form = ~ 1 | Cropping.system),
        data = dat2_time,
        na.action = na.omit)
      # check for normality
      resid <- residuals(fit)
      resp_list[[t]]$normality_pvalue[i] <- shapiro.test(resid)$p.value
      resp_list[[t]]$normality_W[i] <- shapiro.test(resid)$statistic
      resp_list[[t]]$transformed[i] <- "yes"
    }
    
    # export model stats
    mod.out <- as.data.frame(anova(fit, type="marginal")) 
    mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
    resp_list[[t]][i,c(4:6)] <- paste0("F=", round(mod.out2$`F-value`,2), ", p=", mod.out2$`p-value`)[2:4]
    resp_list[[t]]$rsq[i] <- paste0(round(r.squaredGLMM(fit)[1],2), ", ", round(r.squaredGLMM(fit)[2],2)) # marginal, conditional
    # replace p's with asterisks to denote significance
    ps <- mod.out2$`p-value`[2:4]
    for (k in 1:3){
      if(ps[k] >= 0.05){resp_list[[t]][i, k+8] <- "N.S."}
      else if(ps[k] >= 0.01){resp_list[[t]][i, k+8] <- "*"}
      else if(ps[k] >= 0.001){resp_list[[t]][i, k+8] <- "**"}
      else {resp_list[[t]][i, k+8] <- "***"}
    }

    
    ### tukey: Cover * Cropping system
    test1 <- emmeans(fit, ~ Cover*Cropping.system)
    testlet <- cld(test1, type = "response", Letters = "ABCDEFGH", reversed = TRUE)
    testlet <- testlet %>% dplyr::mutate_each_(funs(factor(.)),c("Cover","Cropping.system"))
    testlet$Cropping.system <- factor(testlet$Cropping.system, levels(testlet$Cropping.system)[c(1,4,3,2)])
    testlet$Cover <- factor(testlet$Cover, levels(testlet$Cover)[c(2,4,1,5,3)])
    testlet <- testlet[order(testlet$Cropping.system, testlet$Cover),]
    testlet$.group <- gsub(" ", "", testlet$.group)
    if(length(unique(testlet$.group))==1) {testlet$.group <- rep(" ", length(testlet$.group))}
    write.csv(testlet, paste0("Model-output/anova_", region, "/cover by cropsys/",response_var, "_", timepoint,".csv"))
    
    ### tukey: Cropping.system 
    test1 <- emmeans(fit, ~ Cropping.system)
    testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
    testlet <- testlet %>% dplyr::mutate_each_(funs(factor(.)),c("Cropping.system"))
    testlet$Cropping.system <- factor(testlet$Cropping.system, levels(testlet$Cropping.system)[c(1,4,3,2)])
    testlet <- testlet[order(testlet$Cropping.system),]
    testlet$.group <- gsub(" ", "", testlet$.group)
    if(length(unique(testlet$.group))==1) {testlet$.group <- rep(" ", length(testlet$.group))}
    write.csv(testlet, paste0("Model-output/anova_", region, "/cropsys/",response_var, "_", timepoint,".csv"))
    
    ### tukey: Cover 
    test1 <- emmeans(fit, ~ Cover)
    testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
    testlet <- testlet %>% dplyr::mutate_each_(funs(factor(.)),c("Cover"))
    testlet$Cover <- factor(testlet$Cover, levels(testlet$Cover)[c(2,4,1,5,3)])
    testlet <- testlet[order(testlet$Cover),]
    testlet$.group <- gsub(" ", "", testlet$.group)
    if(length(unique(testlet$.group))==1) {testlet$.group <- rep(" ", length(testlet$.group))}
    write.csv(testlet, paste0("Model-output/anova_", region, "/cover/",response_var, "_", timepoint,".csv"))
    
  }
}

resp_list[[2]] # check results

# now combine all the output into one dataframe
# When you supply a column name with the `.id` argument, a new
# column is created to link each row to its original data frame
resp <- bind_rows(resp_list, .id = "timepoint")

View(resp) # check
resp <- resp[order(resp$resp),]
write.csv(resp, paste0("Model-output/anova_", region, "/all-anova.csv"))





















