library(gtools)
library(dunn.test)

###################################################################
## Load the data
phi_all <- read.csv("Genetic_data/Phi.fst.filtered.csv", header=T)
ecoregions <- c("CoolTemperate" ,"SouthWest", "WarmTemperate","SubTropical")

### Build the ecoregion table
results <- as.data.frame(matrix(NA, 4,4, 
                                dimnames=list(ecoregions,
                                              c("sp_all", "sp_sign", "fst_all", "fst_sign"))))

for (e in 1:4) {
  phi_e <- phi_all[which(phi_all$Bioregion == ecoregions[e]),]
  phi_e_sign <- phi_e[which(phi_e$Pval<0.05),]
  
  results[e,1]<-nrow(phi_e)
  results[e,2]<-nrow(phi_e_sign)
  results[e,3]<-mean(phi_e$Phi_fst)
  results[e,4]<-mean(phi_e_sign$Phi_fst)
}

write.csv(results, file="Genetic_data/Phi_fst_per_bioregion.csv")

#########################
#Plot the results
#########################
tiff(file="Figure4_phifst.tiff", width = 12, height = 25, units = "cm", res=300)
layout(matrix(1:2,2,1,F)) ; layout.show(2)

par(mar=c(2, 5, 3, 1))  # c(bottom, left, top, right)

barplot( results$sp_all, ylim=c(0,6), col="black", beside=FALSE, axes=F) 
par(new=TRUE) 
barplot( results$sp_sign, ylim=c(0,6), col="gray50", beside=FALSE,
         cex.axis=1.6, cex.lab=1.6, cex.names=1.5,
         ylab="Number of species") 


par(mar=c(4, 5, 1, 1))  # c(bottom, left, top, right)

barplot( results$fst_sign, ylim=c(0,0.3), col="gray50", beside=FALSE, axes=F) 
par(new=TRUE) 
barplot( results$fst_all, ylim=c(0,0.3), col="black", beside=FALSE,
         cex.axis=1.6, cex.lab=1.6, cex.names=1.5,
         names.arg=c("CT", "SW", "WT", "ST"),
         ylab=expression('Mean '*phi['ST']),
         xlab="Bioregions") 

legend("topright", 
       legend = c("All", expression('Significant '*phi['ST'])), 
       fill = c("black", "gray50"), ncol = 1, cex=1.3)


dev.off()

##########################
##########################
# Krusakll wallis test for difference between ecoregions
phi_sign <- phi_all[which(phi_all$Pval<0.05),]

kw_all<-kruskal.test(phi_all$Phi_fst ~ factor(phi_all$Bioregion)) # not signif
kw_sign<-kruskal.test(phi_sign$Phi_fst ~ phi_sign$Bioregion) # not signif

dunn.test(phi_all$Phi_fst, factor(phi_all$Bioregion),  method="bonferroni", 
          table=F,list=T, alpha=0.05) 
dunn.test(phi_sign$Phi_fst, factor(phi_sign$Bioregion))
