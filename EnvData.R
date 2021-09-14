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
layersTerr <- list_layers("WorldClim")

# print the Bio-ORACLE citation
print(dataset_citations("Bio-ORACLE"))

# print the citation for a MARSPEC paleo layer
print(layer_citations("MS_biogeo02_aspect_NS_21kya"))
print(layer_citations("WC_bio1_cclgm"))


######################## DOWNLOAD #############################
######## CONTEMPORARY

# download SST mean & range
# Precipitations mean (bio12), wettest month (bio13) and driest month (bio14)
load_layers( layercodes = c("BO2_tempmean_ss", "BO2_temprange_ss",
                            "WC_bio12", "WC_bio13", "WC_bio14") , 
             equalarea=FALSE, rasterstack=TRUE,
             datadir="Env_Layers_contemporary")



######## HISTORICAL
# List of paleoterrestrial recors
paleo.ter<-list_layers_paleo(datasets = c("WorldClim"), epoch = "Last Glacial Maximum")
#can choose cc or mr models
load_layers( layercodes = c("WC_bio14_cclgm", "WC_bio12_cclgm", "WC_bio13_cclgm") ,
             equalarea=FALSE, rasterstack=TRUE,
             datadir="Env_Layers_historical")

## list the available paleo marine recors
paleo.mar<-list_layers_paleo(datasets = c("MARSPEC"), epoch = "Last Glacial Maximum")
load_layers( layercodes = c( "MS_biogeo13_sst_mean_21kya_noCCSM",
                             "MS_biogeo16_sst_range_21kya_noCCSM") , 
             equalarea=FALSE, rasterstack=TRUE,
             datadir="Env_Layers_historical")


######################## WORK WITH LAYERS #############################
# Load marine contemporary rasters from the directory
contemp.mar<- load_layers( layercodes = c("BO2_tempmean_ss", "BO2_temprange_ss") , 
                           equalarea=FALSE, rasterstack=TRUE,
                           datadir="Env_Layers_contemporary")


# Crop raster to fit South Africa 
sa.ext <- extent(14, 35, -37, 21) 
contemp.mar.crop <- crop(contemp.mar, sa.ext) 

# Generate a nice color ramp and plot the map 
my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127")) 
plot(contemp.mar.crop[[1]],col=my.colors(1000),axes=FALSE, box=FALSE) 
title(cex.sub = 1.25, sub = "Mean SST (ºC)") 

############################# EXTRACT DATA FOR SITES ######################################
# Load the data
contemp.mar<- load_layers( layercodes = c("BO2_tempmean_ss", "BO2_temprange_ss") , 
                           equalarea=FALSE, rasterstack=TRUE,
                           datadir="Env_Layers_contemporary")

contemp.terr <- load_layers( layercodes = c("WC_bio12","WC_bio13", "WC_bio14") , 
                             equalarea=FALSE, rasterstack=TRUE,
                             datadir="Env_Layers_contemporary")

paleo.terr <-load_layers( layercodes = c("WC_bio12_cclgm", 
                                         "WC_bio13_cclgm",
                                         "WC_bio14_cclgm"),
                          equalarea=FALSE, rasterstack=TRUE,
                          datadir="Env_Layers_historical")

paleo.mar <- load_layers( layercodes = c( "MS_biogeo13_sst_mean_21kya_noCCSM",
                                          "MS_biogeo16_sst_range_21kya_noCCSM")  , 
                          equalarea=FALSE, rasterstack=TRUE,
                          datadir="Env_Layers_historical")

# load population coordinates
h_df <- read.csv("Genetic_data/hap.div.filtered.csv",  header = T)
n_df <- read.csv("Genetic_data/nuc.div.filtered.csv", header = T)


# Build result matrices (1 per dataset)
env1<-env3<- matrix(NA,nrow(n_df), 2)
env2<-env4<- matrix(NA,nrow(n_df), 3)

# extract values of the sampling sites for each raster
for (i in 1:2) {
  env1[,i] <- extract(contemp.mar[[i]], n_df[,c("Long", "Lat")], 
                      method='bilinear', df=T)[,2]
}

for (i in 1:3) {
  env2[,i] <- extract(contemp.terr[[i]], n_df[,c("Long", "Lat")], 
                      method='bilinear', df=T)[,2]
}

for (i in 1:2) {
  env3[,i] <- extract(paleo.mar[[i]], n_df[,c("Long", "Lat")], 
                      buffer=100000, fun=mean,df=T)[,2]
}

for (i in 1:3) {
  env4[,i] <- extract(paleo.terr[[i]], n_df[,c("Long", "Lat")], 
                      method='bilinear', df=T)[,2]
}

## Create the precipitation range by making the difference between min and max
env2 <- as.data.frame(env2) %>%
  mutate(PrecRange = V2 - V3)

env4 <- as.data.frame(env4) %>%
  mutate(PrecRange = V2 - V3)

## Combine the data in one dataframe
my.sites.environment <- cbind.data.frame(Name=n_df[, c("Sites", "Lat", "Long")] ,env1,env2[,c("V1", "PrecRange")],env3,env4[,c("V1", "PrecRange")]) 
colnames(my.sites.environment) <- c("Sites","Lat","Long",
                                    "SSTmean_cont", "SSTrange_cont",
                                    "PrecMean_cont", "PrecRange_cont",
                                    "SSTmean_paleo", "SSTrange_paleo",
                                    "PrecMean_paleo", "PrecRange_paleo") 

write.csv(my.sites.environment, file="Data_Enviro_Sites_nucl.div_20210914.csv", row.names=F)


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



