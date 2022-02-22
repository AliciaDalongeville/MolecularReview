#########################################################################
## Scripts from Bio-Oracle to extract environmental Data
# http://www.bio-oracle.org/code.php

## sdmpredictors tutorial: 
##https://cran.r-project.org/web/packages/sdmpredictors/vignettes/quickstart.html
##########################################################################
#install.packages("sdmpredictors")
#install.packages("leaflet")

# Load package 
library(sdmpredictors) 
library(dplyr)

######################## EXPLORE #############################
# exploring the marine datasets 
datasets <- list_datasets(terrestrial = F, marine = T)
# MARSPEC Dataset countains historical data

# List the  layers 
layers <- list_layers(datasets)

# print the Bio-ORACLE citation
print(dataset_citations("Bio-ORACLE"))

# print the citation for a MARSPEC paleo layer
print(layer_citations("MS_biogeo02_aspect_NS_21kya"))
print(layer_citations("WC_bio1_cclgm"))


######################## DOWNLOAD #############################
######## CONTEMPORARY

# download SST mean & range, SSS mean & range
load_layers( layercodes = c("BO2_tempmean_ss", "BO2_temprange_ss",
                            "BO2_salinitymean_ss", "BO2_salinityrange_ss") , 
             equalarea=FALSE, rasterstack=TRUE,
             datadir="Env_Layers_contemporary")



######## HISTORICAL
## list the available paleo marine recors
paleo.mar<-list_layers_paleo(datasets = c("MARSPEC"), epoch = "mid-Holocene")
load_layers( layercodes = c( "MS_biogeo13_sst_mean_21kya_noCCSM",
                             "MS_biogeo16_sst_range_21kya_noCCSM",
                             "MS_biogeo08_sss_mean_21kya_noCCSM",
                             "MS_biogeo11_sss_range_21kya_noCCSM",
                             "MS_biogeo13_sst_mean_6kya",
                             "MS_biogeo16_sst_range_6kya",
                             "MS_biogeo08_sss_mean_6kya",
                             "MS_biogeo11_sss_range_6kya" ) , 
             equalarea=FALSE, rasterstack=TRUE,
             datadir="Env_Layers_historical")


######################## WORK WITH LAYERS #############################
# Load marine contemporary rasters from the directory
contemp.mar<- load_layers( layercodes = c("BO2_tempmean_ss", "BO2_temprange_ss",
                                          "BO2_salinitymean_ss", "BO2_salinityrange_ss") , 
                           equalarea=FALSE, rasterstack=TRUE,
                           datadir="Env_Layers_contemporary")


# Crop raster to fit South Africa 
sa.ext <- extent(14, 35, -37, -21) 
paleo.LGM.crop <- crop(paleo.LGM, sa.ext) 

# Generate a nice color ramp and plot the map 
my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127")) 
plot(paleo.LGM.crop[[1]],col=my.colors(1000),axes=FALSE, box=FALSE) 
title(cex.sub = 1.25, sub = "Mean SST LGM (ºC)") 

############################# EXTRACT DATA FOR SITES ######################################
# Load the data
all_layers <- load_layers( layercodes = c("BO2_tempmean_ss", "BO2_temprange_ss",
                                          "BO2_salinitymean_ss", "BO2_salinityrange_ss",
                                          "MS_biogeo13_sst_mean_21kya_noCCSM",
                                          "MS_biogeo16_sst_range_21kya_noCCSM",
                                          "MS_biogeo08_sss_mean_21kya_noCCSM",
                                          "MS_biogeo11_sss_range_21kya_noCCSM",
                                          "MS_biogeo13_sst_mean_6kya",
                                          "MS_biogeo16_sst_range_6kya",
                                          "MS_biogeo08_sss_mean_6kya",
                                          "MS_biogeo11_sss_range_6kya") , 
                           equalarea=FALSE, rasterstack=TRUE,
                           datadir="Env_Layers_contemporary")

# LGM layers
paleo.LGM <- load_layers( layercodes = c( "MS_biogeo13_sst_mean_21kya_noCCSM",
                                          "MS_biogeo16_sst_range_21kya_noCCSM",
                                          "MS_biogeo08_sss_mean_21kya_noCCSM",
                                          "MS_biogeo11_sss_range_21kya_noCCSM")  , 
                          equalarea=FALSE, rasterstack=TRUE,
                          datadir="Env_Layers_historical/LGM")

# mid-Holocene layers
paleo.MH <- load_layers( layercodes = c( "MS_biogeo13_sst_mean_6kya",
                                         "MS_biogeo16_sst_range_6kya",
                                         "MS_biogeo08_sss_mean_6kya",
                                         "MS_biogeo11_sss_range_6kya")  , 
                          equalarea=FALSE, rasterstack=TRUE,
                          datadir="Env_Layers_historical/MH")

# load population coordinates
h_df <- read.csv("Genetic_data/hap.div.filteredNEW.csv",  header = T)
n_df <- read.csv("Genetic_data/nuc.div.filteredNEW.csv", header = T)


# Build result matrices (1 per dataset)
env1<-env2<-env3 <- matrix(NA,nrow(h_df), 4)

# extract values of the sampling sites for each raster
for (i in 1:4) {
  env1[,i] <- raster::extract(contemp.mar[[i]], h_df[,c("Long", "Lat")], 
                      method='bilinear', df=T)[,2]
  env2[,i] <- raster::extract(paleo.LGM[[i]], h_df[,c("Long", "Lat")], 
                      method='bilinear', df=T)[,2]
  env3[,i] <- raster::extract(paleo.MH[[i]], h_df[,c("Long", "Lat")], 
                      method='bilinear', df=T)[,2]
}


## Combine the data in one dataframe
my.sites.environment <- cbind.data.frame(Name=h_df[, c("Site", "Lat", "Long")] ,env1,env2,env3) 
colnames(my.sites.environment) <- c("Site","Lat","Long",
                                    "SSTmean_cont", "SSTrange_cont",
                                    "SSSmean_cont", "SSSrange_cont",
                                    "SSTmean_LGM", "SSTrange_LGM",
                                    "SSSmean_LGM", "SSSrange_LGM",
                                    "SSTmean_MH", "SSTrange_MH",
                                    "SSSmean_MH", "SSSrange_MH") 

write.csv(my.sites.environment, file="Data_Enviro_Sites_hap.div_20210925.csv", row.names=F)


######################################
# Extract mean values per bioregion for the Phi_fst analysis
library(raster)
library(maptools)
library(rgdal)

regions = rgdal::readOGR("/Users/alicia/Documents/Stellenbosch_post-doc/GIS Layers/SA_bioregions_enviro.shp")

env1<- env3 <- matrix(NA,5, 2)
env2<-env4<- matrix(NA,5, 3)


for (i in 1:2) {
  v <- extract(contemp.mar[[i]], regions)
  # mean for each polygon
  env1[,i] <- unlist(lapply(v, function(x) if (!is.null(x)) mean(x, na.rm=TRUE) else NA ))
}

for (i in 1:3) {
  v <- extract(contemp.terr[[i]], regions)
  env2[,i] <- unlist(lapply(v, function(x) if (!is.null(x)) mean(x, na.rm=TRUE) else NA ))
}

for (i in 1:2) {
  v <- extract(paleo.mar[[i]], regions)
  env3[,i] <- unlist(lapply(v, function(x) if (!is.null(x)) mean(x, na.rm=TRUE) else NA ))
}

for (i in 1:3) {
  v <- extract(paleo.terr[[i]], regions)
  env4[,i] <- unlist(lapply(v, function(x) if (!is.null(x)) mean(x, na.rm=TRUE) else NA ))
}

## Create the precipitation range by making the difference between min and max
env2 <- as.data.frame(env2) %>%
  mutate(PrecRange = V2 - V3)

env4 <- as.data.frame(env4) %>%
  mutate(PrecRange = V2 - V3)

## Combine the data in one dataframe
bioreg.environment <- cbind.data.frame(regions$Bioregion,env1,env2[,c("V1", "PrecRange")],env3,env4[,c("V1", "PrecRange")]) 

colnames(bioreg.environment) <- c("Bioregions","SSTmean_cont", "SSTrange_cont",
                                  "PrecMean_cont", "PrecRange_cont",
                                  "SSTmean_paleo", "SSTrange_paleo",
                                  "PrecMean_paleo", "PrecRange_paleo") 

write.csv(bioreg.environment, file="Data_Enviro_Bioregions_20210914.csv", row.names=F)

################
## Compute climate stability (difference between SST/SSS between present and LGM/MH)
env_n <- read.csv("Data_Enviro_Sites_nucl.div_20211005.csv", header=T)
env_h <- read.csv("Data_Enviro_Sites_hap.div_20211005.csv", header=T)

## Calculate stability
env_n <- env_n %>%
  mutate(SST_stability_LGM =  abs(SSTmean_cont - SSTmean_LGM), .after = SSTrange_LGM) %>%
  mutate(SSS_stability_LGM =  abs(SSSmean_cont - SSSmean_LGM), .after = SSSrange_LGM) %>%
  mutate(SST_stability_MH =  abs(SSTmean_cont - SSTmean_MH), .after = SSTrange_MH) %>%
  mutate(SSS_stability_MH =  abs(SSSmean_cont - SSSmean_MH), .after = SSSrange_MH)

env_h <- env_h %>%
  mutate(SST_stability_LGM =  abs(SSTmean_cont - SSTmean_LGM), .after = SSTrange_LGM) %>%
  mutate(SSS_stability_LGM =  abs(SSSmean_cont - SSSmean_LGM), .after = SSSrange_LGM) %>%
  mutate(SST_stability_MH =  abs(SSTmean_cont - SSTmean_MH), .after = SSTrange_MH) %>%
  mutate(SSS_stability_MH =  abs(SSSmean_cont - SSSmean_MH), .after = SSSrange_MH)

## Write results
write.csv(env_h, file="Data_Enviro_Sites_hap.div_20211005.csv", row.names=F)
write.csv(env_n, file="Data_Enviro_Sites_nucl.div_20211005.csv", row.names=F)

## Calculate pairwise correlations between predictors 
M <- cor(env_h[,4:19], use="pairwise.complete.obs")

# Plot the correlogram
png("Env_predictors_correlogram.png", width = 800, height = 800)
corrplot::corrplot(M, type="lower", diag=F,
                   method = "color", rect.col = "black", 
                   rect.lwd = 5,outline = F,cl.cex = 0.25, 
                   number.cex = 0.8,addCoef.col = "black",tl.col = "black",
                   col = colorRampPalette(c("red", "white","midnightblue"))(200))
dev.off()
