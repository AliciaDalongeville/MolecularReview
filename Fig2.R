library(gtools)
library(FSA)
library(dunn.test)
library(Hmisc)
library(dplyr)

# load population coordinates
h_df <- read.csv("Genetic_data/hap.div.filteredNEW.csv",  header = T)
n_df <- read.csv("Genetic_data/nuc.div.filteredNEW.csv", header = T)

###########################################################
## Linear regression
X_h <- h_df[,c("Bioregion","Substrate")]
X_n <- n_df[,c("Bioregion","Substrate")]

Y_h<-logit(h_df$h) # logit to normalise
Y_n<-logit(n_df$pi)

data_h <- cbind(Y_h,X_h)
colnames(data_h)[1] <- "h"
data_n <- cbind(Y_n,X_n)
colnames(data_n)[1] <- "pi"

## Kruskal-Wallis Test to compare haplotype diversity per traits
res_kw <- matrix(NA,4,2, dimnames=list(c("h_Bioregion", "h_Habitat", "pi_Bioregion", "pi_Habitat"), 
                                       c("chi-sq", "pval")))

# h for Bioregions
kw_h_bio <- kruskal.test(Y_h ~ X_h[,"Bioregion"]) 
res_kw[1,] <- c(kw_h_bio$statistic,kw_h_bio$p.value)
# h for Habitat
kw_h_hab <- kruskal.test(Y_h ~ X_h[,"Substrate"]) 
res_kw[2,] <- c(kw_h_hab$statistic,kw_h_hab$p.value)
# pi for Bioregions
kw_n_bio <- kruskal.test(Y_n ~ X_n[,"Bioregion"]) 
res_kw[3,] <- c(kw_n_bio$statistic,kw_n_bio$p.value)
# pi for Habitat
kw_n_hab <- kruskal.test(Y_n ~ X_n[,"Substrate"]) 
res_kw[4,] <- c(kw_n_hab$statistic,kw_n_hab$p.value)

signif <- res_kw[which(res_kw[,2] < 0.05),] # Only bioregion is significant
write.csv(res_kw, file="Resultats_KWtest.csv")

# Dunn posthoc test to test significance of pairwise differences between bioregions
dunn_h = dunn.test(Y_h,  X_h[,"Bioregion"],  method="bonferroni", 
                   table=F,list=T, alpha=0.1) 
res_dunn_h <- cbind.data.frame(dunn_h$comparisons,dunn_h$Z, dunn_h$P.adjusted)
dunn.signif_h <- res_dunn_h[which(res_dunn_h[,3] < 0.05),]
write.csv(res_dunn_h, file="Resultats_Dunn_h_Bioregions.csv")

dunn_n = dunn.test(Y_n,  X_n[,"Bioregion"],  method="bonferroni", 
                   table=F,list=T, alpha=0.1) 
res_dunn_n <- cbind.data.frame(dunn_n$comparisons,dunn_n$Z, dunn_n$P.adjusted)
dunn.signif_n <- res_dunn_n[which(res_dunn_n[,3] < 0.05),]
write.csv(res_dunn_n, file="Resultats_Dunn_pi_Bioregions.csv")

#### Plot the dots chart
## See tutorial http://www.sthda.com/english/wiki/dot-charts-r-base-graphs

Y_h<-h_df$h
Y_n<-n_df$pi

data_h <- cbind(Y_h,X_h)
colnames(data_h)[1] <- "h"
data_n <- cbind(Y_n,X_n)
colnames(data_n)[1] <- "pi"

## Calculate the mean and standard errors per group
plot_df <- as.data.frame(matrix(NA,12,4))
nb_mod<- c(4,2)
plot_df[,1]<- c(rep(colnames(X_h), nb_mod), rep(colnames(X_n), nb_mod))
e=1

# define standard error function
st.err <- function(x) {
  sd(x)/sqrt(length(x))
}

# calculate mean and se for h
# Bioregions
plot_df[1:4,2] <- as.character(aggregate(Y_h ~ X_h[,1] , data_h,  mean)[,1])
plot_df[1:4,3]<-aggregate(Y_h ~ X_h[,1] , data_h,  mean)[,2]
plot_df[1:4,4]<-aggregate(Y_h ~ X_h[,1] , data_h,  st.err)[,2]
# Habitat
plot_df[5:6,2] <- as.character(aggregate(Y_h ~ X_h[,2] , data_h,  mean)[,1])
plot_df[5:6,3]<-aggregate(Y_h ~ X_h[,2] , data_h,  mean)[,2]
plot_df[5:6,4]<-aggregate(Y_h ~ X_h[,2] , data_h,  st.err)[,2]

# calculate mean and se for pi
# Bioregions
plot_df[7:10,2] <- as.character(aggregate(Y_n ~ X_n[,1] , data_n,  mean)[,1])
plot_df[7:10,3]<-aggregate(Y_n ~ X_n[,1] , data_n,  mean)[,2]
plot_df[7:10,4]<-aggregate(Y_n ~ X_n[,1] , data_n,  st.err)[,2]
# Habitat
plot_df[11:12,2] <- as.character(aggregate(Y_n ~ X_n[,2] , data_n,  mean)[,1])
plot_df[11:12,3]<-aggregate(Y_n ~ X_n[,2] , data_n,  mean)[,2]
plot_df[11:12,4]<-aggregate(Y_n ~ X_n[,2] , data_n,  st.err)[,2]

colnames(plot_df)<- c("Variable", "Mod", "mean", "st.err")
write.csv(plot_df, file="Mean_sd_h_pi_per_variable.csv")

##########################
# PLOT
###########################
plot_df_ALL <- read.csv("Mean_sd_h_pi_per_variable.csv", header=T)

## create pch vector: different symbols when dunn signif
pch<-c(21,21,21,24,21,21)

fill<-c("white", "black", "black", "black",  "white", "white")

labels <- c("Cool Temperate" , "South West" , "Warm Temperate", "Sub-Tropical", 
            "Hard substrate","Soft substrate")

## Draw plot
SE_h=plot_df_ALL$st.err[1:6]
SE_pi=plot_df_ALL$st.err[7:12]

plot_h <- plot_df_ALL[1:6,]
plot_pi <- plot_df_ALL[7:12,]


tiff(file="Fig2_a.tiff", width = 28, height = 18, units = "cm", res=300)
layout(matrix(1:2,1,2,F)) ; layout.show(2)

# Plot h
par(mar=c(3, 9, 1, 0))  # c(bottom, left, top, right)
dotchart2(plot_h$mean, labels = labels, groups = plot_h$Variable,
          xlab="", cex.labels=1.2, xlim=c(0.55,1), cex=1.2,
          cex.group.labels=1.2, groupfont=2, leavepar=TRUE, lty=3, 
          lcolor="gray60", pch=pch, col="black", bg=fill, 
          dotsize=1.3, sort.=F)

dotchart2(plot_h$mean+SE_h, labels = labels, groups = plot_h$Variable, 
          pch="|", dotsize=0.9, add=T,leavepar=TRUE, reset.par=F, sort.=F) 
dotchart2(plot_h$mean-SE_h, labels = labels, groups = plot_h$Variable, 
          pch="|",dotsize=0.9,  add=T,leavepar=TRUE, reset.par=F,sort.=F)

mtext(side = 1, line = 2, 'h', cex=1.5)
dev.off()

# Plot pi
tiff(file="Fig2_b.tiff", width = 14, height = 18, units = "cm", res=300)
par(mar=c(3, 2, 1, 7))  # c(bottom, left, top, right)

dotchart2(plot_pi$mean, labels = "", groups = plot_pi$Variable,
          xlab='', ylab='', cex.labels=1.2, xlim=c(0,0.02), cex=1.2,
          cex.group.labels=1.2, groupfont=0, leavepar=TRUE, lty=3, 
          lcolor="gray60", pch=pch, col="black", bg=fill, 
          dotsize=1.3, sort.=F)

dotchart2(plot_pi$mean+SE_pi, labels = labels, groups = plot_pi$Variable, 
          pch="|", dotsize=0.9, add=T,leavepar=TRUE, reset.par=F, sort.=F) 
dotchart2(plot_pi$mean-SE_pi, labels = labels, groups = plot_pi$Variable, 
          pch="|",dotsize=0.9,  add=T,leavepar=TRUE, reset.par=F,sort.=F)

mtext(side = 1, line = 2, expression(pi), cex=1.5)

dev.off()

#############################

#######################################################################################
## Heatmap of hap and nuc diversity
#######################################################################################
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
library(ggplot2)
library(mapdata)
library(dplyr)
library(viridis)
library(colorBlindness)
library(spData)
library(rgeos)
library(RColorBrewer)
library(ggthemes)
library(gridExtra)


# load population coordinates
h_df <- read.csv("Genetic_data/hap.div.filteredNEW.csv",  header = T)
n_df <- read.csv("Genetic_data/nuc.div.filteredNEW.csv", header = T)

#Import present day coastline shapefile
my_spdf <- readOGR("Env_Layers_historical/LGM/Extrapolation/AFbuff2.shp")
### may want to clip this down as it extends over much of central Africa
sa.ext <- extent(14, 35, -37, -27) 
my_spdf <- raster::crop(my_spdf, sa.ext)

#Transform shapefile into empy raster
pgeo <- spTransform(my_spdf, CRS('+proj=longlat +datum=WGS84 +no_defs'))
ext <- floor(extent(pgeo))
rr <- raster(ext, res=0.1)
rr <- rasterize(pgeo, rr, field=1)
crs(rr) <- '+proj=longlat +datum=WGS84 +no_defs'
plot(rr)

xy_h <- h_df[,c("Lat", "Long")]
xy_n <- n_df[,c("Lat", "Long")]
d_h <- h_df[,c("Lat", "Long", "h")]
d_n <- n_df[,c("Lat", "Long", "pi")]

dsp_h <- SpatialPoints(xy_h, proj4string=CRS('+proj=longlat +datum=WGS84 +no_defs'))
dsp_n <- SpatialPoints(xy_n, proj4string=CRS('+proj=longlat +datum=WGS84 +no_defs'))
dsp_h <- SpatialPointsDataFrame(dsp_h, d_h)
dsp_n <- SpatialPointsDataFrame(dsp_n, d_n)
dta_h <- spTransform(dsp_h, CRS('+proj=longlat +datum=WGS84 +no_defs'))
dta_n <- spTransform(dsp_n, CRS('+proj=longlat +datum=WGS84 +no_defs'))


#± Fit the model for interpolation 
#gs_h <- gstat(formula=h~1, locations=dta_h, nmax=10, nmin=3, set = list(idp = 0.5))

# Thin Plate Spline Regression
fit_TPS_h <- fields::Tps( # using {fields}
  x = as.matrix(d_h[, c("Long", "Lat")]), # accepts points but expects them as matrix
  Y = d_h$h,  # the dependent variable
  miles = FALSE     # EPSG 25833 is based in meters
)
fit_TPS_n <- fields::Tps( # using {fields}
  x = as.matrix(d_n[, c("Long", "Lat")]), # accepts points but expects them as matrix
  Y = d_n$pi,  # the dependent variable
  miles = FALSE     # EPSG 25833 is based in meters
)

# Perform interpolation
idw_h <- interpolate(rr,fit_TPS_h)
idw_n <- interpolate(rr,fit_TPS_n)
## [inverse distance weighted interpolation]
idwmsk_h <- mask(idw_h, rr)
idwmsk_n <- mask(idw_n, rr)


###################
# Generate  the map 
wH <- map_data("worldHires",  xlim=c(14,35), ylim=c(-37,-21)) # subset polygons surrounding SA

# raster as dataframes 
coords_h <- xyFromCell(idwmsk_h, seq_len(ncell(idwmsk_h)))
ndvi_h <- stack(as.data.frame(getValues(idwmsk_h)))
names(ndvi_h) <- c('h', 'variable')
ndvi_h <- cbind(coords_h, ndvi_h)

coords_n <- xyFromCell(idwmsk_n, seq_len(ncell(idwmsk_n)))
ndvi_n <- stack(as.data.frame(getValues(idwmsk_n)))
names(ndvi_n) <- c('pi', 'variable')
ndvi_n <- cbind(coords_n, ndvi_n)


# Change the breaks for the legend 
myfuns <- list(Minimum = min, Maximum = max)
ls_val_h <- unlist(lapply(myfuns, function(f) f(ndvi_h$h, na.rm=T)))
ls_val_n <- unlist(lapply(myfuns, function(f) f(ndvi_n$pi, na.rm=T)))
names(ls_val_h) <- round(ls_val_h,3)
#names(ls_val_n) <- round(ls_val_n,3)
names(ls_val_n) <- c("0.000","0.011")

# Plot haplotype diversity 
map_h <- ggplot(ndvi_h) +  
  geom_polygon(data=wH, aes(x=long, y = lat, group = group), 
               fill="grey80", color="grey50", size=0.25, alpha=0.5) +
  geom_raster(aes(x, y, fill = h)) + 
  scale_fill_viridis(na.value="transparent", breaks=ls_val_h) +
  coord_fixed(xlim=c(14,35), ylim=c(-37,-21), ratio=1.2)+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.title = element_blank())+
  theme(legend.position="bottom") +
  theme(legend.key.width=unit(1.2, "cm")) +
  theme(legend.title = element_blank()) +
  ggtitle("Haplotype diversity") +
  theme(plot.title = element_text(color="black", size=14, face="bold", hjust=0.5)) +
  theme(axis.text.x=element_text(colour="black",size=12))+
  theme(axis.text.y=element_text(colour="black", size=12))+
  theme(legend.text=element_text(colour="black",size=12))

# Plot nucleotide diversity
map_n <- ggplot(ndvi_n) +  
  geom_polygon(data=wH, aes(x=long, y = lat, group = group), 
               fill="grey80", color="grey50", size=0.25, alpha=0.5) +
  geom_raster(aes(x, y, fill = pi)) + 
  scale_fill_viridis(na.value="transparent", breaks=ls_val_n) +
  coord_fixed(xlim=c(14,35), ylim=c(-37,-21), ratio=1.2)+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.title = element_blank())+
  theme(legend.position="bottom") +
  theme(legend.key.width=unit(1.2, "cm")) +
  theme(legend.title = element_blank()) +
  ggtitle("Nucleotide diversity") +
  theme(plot.title = element_text(color="black", size=14, face="bold", hjust=0.5)) +
  theme(axis.text.x=element_text(colour="black",size=12))+
  theme(axis.text.y=element_text(colour="black", size=12))+
  theme(legend.text=element_text(colour="black",size=12))

tiff(file="Heatmaps_diversity.tiff", width = 28, height = 18, units = "cm", res=300)
grid.arrange(map_h,map_n,nrow=1,ncol=2) #use package gridExtra   
dev.off()

##############################################################################
## Data for Table 2
h_df %>%
  group_by(Bioregion,Substrate) %>%
  summarise(mean=mean(h), sd=sd(h),n_sites=n(),n_species=n_distinct(spp))

h_df %>%
  group_by(Bioregion) %>%
  summarise(mean=mean(h), sd=sd(h),n_sites=n(),n_species=n_distinct(spp))

h_df %>%
  group_by(Substrate) %>%
  summarise(mean=mean(h), sd=sd(h),n_sites=n(),n_species=n_distinct(spp))

n_df %>%
  group_by(Bioregion,Substrate) %>%
  summarise(mean=mean(pi), sd=sd(pi))

n_df %>%
  group_by(Bioregion) %>%
  summarise(mean=mean(pi), sd=sd(pi))

n_df %>%
  group_by(Substrate) %>%
  summarise(mean=mean(pi), sd=sd(pi))
