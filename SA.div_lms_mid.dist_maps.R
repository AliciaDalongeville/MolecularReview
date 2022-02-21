## SA div paper - linear models, plotting sea surface temps, & calculate distance from core ###

library(lme4)
library(plyr)
library("dplyr")
library("PerformanceAnalytics")
library(glmmTMB)
library(r2glmm)
library(glmm)


normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

sa_div_data_all_norm$h_norm<-normalize(sa_div_data_all_norm$h)
sa_div_data_all_norm$pi_norm<-normalize(sa_div_data_all_norm$pi)

library(lmerTest)
m1.pi <- lmer(pi_norm ~ Range_position + SSTmean_cont + SSSmean_cont + SSTrange_cont + SST_stability_MH + (1 | spp),  data = sa_nuc_div_norm_out_filtered)
            
m1.h <- lmer(h_norm ~ Range_position + SSTmean_cont + SSSmean_cont + SSTrange_cont + SST_stability_MH + (1 | spp),  data = sa_hap_div_norm_out_filtered)

pdf("hap.div.sst.range.pdf")
sa_hap_div_norm_out_filtered$fit <- predict(m1.h)
ggplot(sa_hap_div_norm_out_filtered,aes(SSTrange_cont, h_norm, col=spp)) + 
    geom_line(aes(y=fit), size=0.8) +
  geom_point(alpha = 0.3) + 
  geom_hline(yintercept=0, linetype="dashed") +
  labs(x = "SSTrange", y = "Haplotype diversity")+
  theme_bw()
dev.off()

pdf("nuc.div.sst.range.pdf")
sa_nuc_div_norm_out_filtered$fit <- predict(m1.pi)
ggplot(sa_nuc_div_norm_out_filtered,aes(SSTrange_cont, pi_norm, col=spp)) + 
  geom_line(aes(y=fit), size=0.8) +
  geom_point(alpha = 0.3) + 
  geom_hline(yintercept=0, linetype="dashed") +
  labs(x = "SSTrange", y = "Nucleotide diversity")+
  theme_bw()
dev.off()

library(sjPlot)

pdf("hap.div.re.pdf")
plot_model(m1.h, type = "re")
dev.off()

pdf("nuc.div.re.pdf")
plot_model(m1.pi, type = "re")
dev.off()

tab_model(m1.h, m1.pi)

#### Code used to create bottom 2 maps in Fig. 1 (top 2 made in QGIS)
library(rgdal)
library(raster)
library(tmap)

setwd("~/Desktop/PhD_stuffies/SDMS/SDM.shps")
map1 <- readOGR("Africa.shp")
SA.sg.ext <- extent(14, 33, -35.5, -25)
map4 <- crop(map1, SA.sg.ext)

sst.range <- 'BO2_temprange_ss_lonlat.tif'
sst.range <- stack(sst.range)
sst.range.c <- crop(sst.range, SA.sg.ext)
sst.r <- map5 + tm_shape(sst.range.c) + tm_raster(style = "cont", palette = "PuRd")+tm_legend(outside=TRUE)

sst.m <- 'BO2_tempmean_ss_lonlat.tif'
sst.m <- stack(sst.m)
sst.m.c <- crop(sst.m, SA.sg.ext)
sst.m <- map5 + tm_shape(sst.m.c) + tm_raster(style = "cont", palette = "YlOrRd")+tm_legend(outside=TRUE)

png("sst.range.mean.png", units="in", width=5, height=5, res=300)
tmap_arrange(sst.r, sst.m, ncol=2)
dev.off()


#### Code trying to calculate dist from center SA genetic diversity paper

library("rgdal")
library("rgeos")

setwd("~/Desktop/PhD_stuffies/SDMS/SDM.shps")
AF.line <- readOGR(dsn="AF.coast.line.shp")
route <- spTransform(AF.line, CRS("+init=epsg:32634"))

wgs84 <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
p = SpatialPointsDataFrame(coords=hap_div_filtered[,5:4], proj4string = wgs84, data=hap_div_filtered)
p = spTransform(p,CRS("+init=epsg:32634"))

d <- sort(gProject(route, p))
diff(d)
#this gives the distance between each pt, have to make sure points are in order to east to west!

#### Script Not Used:

pdf(file="predictor.corrs.pdf") 
chart.Correlation(sa_div_data_all[, c(12:17)], histogram=TRUE, pch=19)
dev.off()

#data_tot=ddply(sa_div_data_all, .(h),mutate,GD_scaled=as.vector(scale(GD))) 
#attach(data_tot)

#data_tot=ddply(sa_div_data_all, c("h", "pi"), transform, x.std = scale(x))

(ggplot(pred.mm) + 
    geom_line(aes(x = x, y = predicted)) +          # slope
    geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                fill = "lightgrey", alpha = 0.5) +  # error band
    geom_point(data = sa_div_data_all_norm,                      # adding the raw data (scaled values)
               aes(x = SSSmean_cont, y = pi_norm, colour = spp)) + 
    labs(x = "Body Length (indexed)", y = "Test Score", 
         title = "Body length does not affect intelligence in dragons") + 
    theme_minimal()
)


sjp.lmer(m0,
         facet.grid = FALSE,
         sort.est = "sort.all")


# Create normal distribution for hap & nuc diversity
#div_data_scale <- sa_div_data_all %>%           
# mutate_at(c("h", "pi"), ~(scale(.) %>% as.vector))


sjp.lmer(m0, vars = "SSTrange_cont", type = "ri.slope")

fixef(m0)

ranef(m0)

#m0 <- glmmTMB( pi_norm ~ Range_position + SSTmean_cont + SSSmean_cont + SSTrange_cont + SST_stability_MH + (1 | spp), family = gaussian,  data=sa_div_data_all)


div_data_cs <- transform(div_data_scale,
RP_cs=scale(Range_position),
SST_m_cs=scale(SSTmean_cont),
SSS_m_cs=scale(SSSmean_cont),
SST_r_cs=scale(SSTrange_cont),
SST_s_cs=scale(SST_stability_MH)
)    