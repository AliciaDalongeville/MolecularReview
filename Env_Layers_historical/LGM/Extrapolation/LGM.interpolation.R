# Script to interpolate LGM raster to present day bathymetry

library(sf)
library(tmap)
library(gstat)
library(readxl)
library(raster)
library(rgdal)
library(fields)
library(automap)
library(dplyr)
library(sdmpredictors) 

#Import present day coastline shapefile
my_spdf <- readOGR("Env_Layers_historical/LGM/Extrapolation/AFbuff2.shp")
### may want to clip this down as it extends over much of central Africa

#Transform shapefile into empy raster
pgeo <- spTransform(my_spdf, CRS('+proj=longlat +datum=WGS84 +no_defs'))
ext <- floor(extent(pgeo))
rr <- raster(ext, res=0.1)
rr <- rasterize(pgeo, rr, field=1)
crs(rr) <- '+proj=longlat +datum=WGS84 +no_defs'
plot(rr)


#Get xy points from LGM raster
paleo.LGM <- load_layers( layercodes = c( "MS_biogeo13_sst_mean_21kya_noCCSM",
                                          "MS_biogeo16_sst_range_21kya_noCCSM",
                                          "MS_biogeo08_sss_mean_21kya_noCCSM",
                                          "MS_biogeo11_sss_range_21kya_noCCSM")  , 
                          equalarea=FALSE, rasterstack=TRUE,
                          datadir="Env_Layers_historical/LGM")

for (l in 1:4) {
  print(names(paleo.LGM)[l])
  layer <- paleo.LGM[[l]]
  # Crop raster to fit South Africa 
  sa.ext <- extent(14, 35, -37, 21) 
  layer.crop <- raster::crop(layer, sa.ext) 
  
  xy <- data.frame(xyFromCell(layer.crop, 1:ncell(layer.crop)))
  v <- getValues(layer.crop)
  i <- !is.na(v)
  xy <- xy[i,]
  v <- v[i]
  d <- cbind(xy,v)
  
  
  dsp <- SpatialPoints(xy, proj4string=CRS('+proj=longlat +datum=WGS84 +no_defs'))
  dsp <- SpatialPointsDataFrame(dsp, d)
  dta <- spTransform(dsp, CRS('+proj=longlat +datum=WGS84 +no_defs'))
  
  #Fit the model for interpolation 
  gs <- gstat(formula=v~1, locations=dta)
  
  
  # Perform interpolation
  idw <- interpolate(rr,gs)
  ## [inverse distance weighted interpolation]
  idwmsk <- mask(idw, rr)
  plot(idwmsk)
  
  # Save the raster
  writeRaster(idwmsk, filename=paste("Env_Layers_historical/LGM/Extrapolation/",names(paleo.LGM)[[l]], "_extrapolated.tif", sep=""), overwrite=T)
}


#####################
## Extract values of hap and nuc div sites
# load population coordinates
h_df <- read.csv("Genetic_data/hap.div.filteredNEW.csv",  header = T)
n_df <- read.csv("Genetic_data/nuc.div.filteredNEW.csv", header = T)

# load the result matrices
res_h <- read.csv("Data_Enviro_Sites_hap.div_20210925.csv", header=T)
res_n <- read.csv("Data_Enviro_Sites_nucl.div_20210925.csv", header=T)

# extract values of the sampling sites for each raster
#first import all files in a single folder as a list 
rastlist <- list.files(path = "Env_Layers_historical/LGM/Extrapolation", pattern='.tif', all.files=TRUE, full.names=TRUE)
# Stack the rasters
paleo.LGM.extra <- stack(rastlist)

index <- c(10,11,8,9)
  for (i in 1:4) {
    e <- index[i]
    res_h[,e] <- raster::extract(paleo.LGM.extra[[i]], h_df[,c("Long", "Lat")], 
                                method='bilinear', df=T)[,2]
    res_n[,e] <- raster::extract(paleo.LGM.extra[[i]], n_df[,c("Long", "Lat")], 
                                     method='bilinear', df=T)[,2]
  }
