#This script runs a redundancy analysis to identify genotype-environment associations, as shown in: Stager M, Senner NR, Swanson DL, Carling MD, Grieves TJ, and Cheviron ZA. 2021. Temperature heterogeneity correlates with intraspecific variation in physiological flexibility in a small endotherm. Nature Communications, accepted.

#This workflow has been modified from Forester et al. 2018's tutorial (available at https://popgen.nescent.org/2018-03-27_RDA_GEA.html) to incorporate variance partitioning

#################################################
library(adegenet)
library(pegas)
library(vegan)
library(factoextra)

SNPS <- read.PLINK("Junco_SNPs.raw") #output from processing RAD data
gen <- as.data.frame(SNPS)
dim(gen)

geo<-read.csv("Source_Data_Fig2.csv") #georeferenced samples listed in Source Data Sheet 2 ("Fig.2")
geo <- geo[geo$Catalog_number!="6056" & geo$Catalog_number!="110508" & geo$Catalog_number!="117690" & geo$Catalog_number!="118252" & geo$Catalog_number!="229176" & geo$Catalog_number!="231996" & geo$Catalog_number!="648103" & geo$Catalog_number!="648671" & geo$Catalog_number!="232004" & geo$Catalog_number!="183387" & geo$Catalog_number!="232018",] #remove 11 individuals that were filtered during RAD processing

sum(is.na(gen)) #total number of sites to impute
sum(is.na(gen)) / (dim(gen)[1]*dim(gen)[2]) #percentage of total

#impute using most common genotype for that population
gen.imp <- apply(gen[geo$Species=="hyemalis",], 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
gen.imp <- rbind(gen.imp, apply(gen[geo$Species=="insularis",], 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x)))))))
gen.imp <- rbind(gen.imp, apply(gen[geo$Species=="phaeonotus",], 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x)))))))
gen.imp <- rbind(gen.imp, apply(gen[geo$Species=="vulcani",], 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x)))))))
gen.imp <- rbind(gen.imp, apply(gen[geo$Species=="bairdi",], 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x)))))))

gen.imp = gen.imp[row.names(gen.imp)[match(row.names(gen), row.names(gen.imp))],]
identical(row.names(gen.imp[row.names(gen.imp)[match(row.names(gen), row.names(gen.imp))],]), geo$Catalog_number) #make sure files are in the same order

sum(is.na(gen.imp)) # No NAs

pca_all <- dudi.pca(gen.imp,scannf=FALSE,scale=FALSE, nf=6)
get_eig(pca_all)[1:6,]

#environmental data
library(raster)

data=getData('worldclim', var='bio', res=2.5, download=TRUE) #bioclim data
r <- data[[c(1:19)]] #using all bioclim variables
names(r) <- c("Mean.Temp", "Diurnal.Range", "Isothermality", "Temp.Seasonality", "Max.Temp", "Min.Temp", "Annual.Range", "Temp.WettestQ", "Temp.DriestQ", "Temp.WarmestQ", "Temp.ColdestQ", "Annual.Precip", "Max.Precip", "Min.Precip", "Precip.Seasonality", "Precip.WettestQ", "Precip.DriestQ", "Precip.WarmestQ", "Precip.ColdestQ") #bioclim variable names
coords <- data.frame(long=geo$Longitude,lat=geo$Latitude)
points <- SpatialPoints(coords, proj4string = r@crs)
values <- extract(r,points) #extract data for sample locales
clim <- cbind.data.frame(coordinates(points),values)
clim $sampleID <- geo$Catalog_number #add catalog number
identical(rownames(gen.imp), clim$sampleID) #in same order as SNP data?
                                
#Plot sampling locations overlaid on Seasonality map (Fig 2a)
count=as.data.frame(table(round(geo$Latitude,1), round(geo$Longitude,1), geo$Plot_color))
count=count[!(count$Freq==0),]

plot(data$bio7, xlim=c(-170,-53), ylim=c(7,75),bty="n", col=gray.colors(10, rev=TRUE)) #Temp Annual Range of NA
points(y=c(as.numeric(as.character(count$Var1))), x=c(as.numeric(as.character(count$Var2))),  pch=21, bg=as.character(count$Var3), col="black", cex=1+(count$Freq/2.5))


library(arm)
clim = data.frame(apply(clim, 2, rescale))

#removing environmental variables that are strongly correlated (r > 0.7)
library(psych)
pairs.panels(clim[,3:21], scale=T, cex=2.5)
pred <- subset(clim, select=-c(long,lat,sampleID,Mean.Temp, Temp.ColdestQ, Min.Temp, Annual.Precip, Temp.Seasonality, Max.Precip, Precip.DriestQ, Isothermality, Precip.WettestQ, Max.Temp, Precip.Seasonality)) #in order of removal, all r < 0.7
pairs.panels(pred, scale=T, cex=2.5)


#run the partial rda
pred$PC1 <- as.vector(pca_all$li[,1])
pred$PC2 <- pca_all$li[,2]

junco.prda <- rda(gen.imp ~ Diurnal.Range + Annual.Range + Temp.WettestQ + Temp.DriestQ + Temp.WarmestQ + Min.Precip + Precip.WarmestQ + Precip.ColdestQ + Condition(PC1+PC2), data=pred)
junco.prda
RsquareAdj(junco.prda) #adjusted R-squared for number of predictors
summary(eigenvals(junco.prda, model = "constrained"))
vif.cca(junco.prda) 

screeplot(junco.prda) #visual representation
signif.full.prda <- anova.cca(junco.prda, parallel=getOption("mc.cores")) # default is permutation=999
signif.full.prda #p = 0.001

signif.axis.prda <- anova.cca(junco.prda, by="axis", parallel=getOption("mc.cores"),  permutations=how(nperm=299)) #may be able to up permutations but intensive
signif.axis.prda #R1 - R3, p = 0.01

#Variance Partitioning
varpart(gen.imp, ~Diurnal.Range, ~ Annual.Range + Temp.WettestQ + Temp.DriestQ + Temp.WarmestQ + Min.Precip + Precip.WarmestQ + Precip.ColdestQ, ~PC1 + PC2, data=pred)

varpart(gen.imp, ~Annual.Range, ~ Diurnal.Range + Temp.WettestQ + Temp.DriestQ + Temp.WarmestQ + Min.Precip + Precip.WarmestQ + Precip.ColdestQ, ~PC1 + PC2, data=pred)

varpart(gen.imp, ~ Temp.WettestQ, ~Diurnal.Range + Annual.Range + Temp.DriestQ + Temp.WarmestQ + Min.Precip + Precip.WarmestQ + Precip.ColdestQ, ~PC1 + PC2, data=pred)

varpart(gen.imp, ~ Temp.DriestQ, ~Diurnal.Range + Annual.Range + Temp.WettestQ + Temp.WarmestQ + Min.Precip + Precip.WarmestQ + Precip.ColdestQ, ~PC1 + PC2, data=pred)

varpart(gen.imp, ~ Temp.WarmestQ, ~Diurnal.Range + Annual.Range + Temp.WettestQ + Temp.DriestQ + Min.Precip + Precip.WarmestQ + Precip.ColdestQ, ~PC1 + PC2, data=pred)

varpart(gen.imp, ~Min.Precip, ~Diurnal.Range + Annual.Range + Temp.WettestQ + Temp.DriestQ + Temp.WarmestQ + Precip.WarmestQ + Precip.ColdestQ, ~PC1 + PC2, data=pred)

varpart(gen.imp, ~Precip.WarmestQ, ~Diurnal.Range + Annual.Range + Temp.WettestQ + Temp.DriestQ+  Temp.WarmestQ + Min.Precip + Precip.ColdestQ, ~PC1 + PC2, data=pred)

varpart(gen.imp, ~Precip.ColdestQ, ~Diurnal.Range + Annual.Range + Temp.WettestQ + Temp.DriestQ+  Temp.WarmestQ + Min.Precip  + Precip.WarmestQ, ~PC1 + PC2, data=pred)

varpart(gen.imp, ~PC1+PC2, ~Diurnal.Range + Annual.Range + Temp.WettestQ + Temp.DriestQ+  Temp.WarmestQ + Min.Precip  + Precip.WarmestQ + Precip.ColdestQ, data=pred)


#Plot PCAs (Fig 2b and Fig2c)
temp <- as.integer(as.factor(paste(geo$Species,geo$Subspecies)))
myCol <- c("gray48","darkolivegreen", "yellowgreen", "darkseagreen1","cadetblue","cyan", "dodgerblue","blue2","purple4","darkorchid1", "plum3", "rosybrown", "coral2", "red", "red4", "darkgoldenrod","darkorange","goldenrod1","yellow2", "bisque", "white")[temp] #alphabetical order

par(mfrow=c(2,3))
plot(pca_all$li[,1:2], col="gray32", bg=myCol, pch=21, cex=2, xlab=paste("PC1 (", round(get_eig(pca_all)[1,2],1), "%)", sep=""), ylab=paste("PC2 (", round(get_eig(pca_all)[2,2],1), "%)", sep=""))

plot(pca_all$li[,c(3,4)], col="gray32", bg=myCol, pch=21, cex=2, xlab=paste("PC3 (", round(get_eig(pca_all)[3,2],1), "%)", sep=""), ylab=paste("PC4(", round(get_eig(pca_all)[4,2],1), "%)", sep=""))

plot.new()
legend("left",fill=unique(myCol)[order(unique(geo$Subspecies))], legend=unique(geo$Subspecies)[order(unique(geo$Subspecies))], bty="n", cex=1.2,text.font=3)

#Plot RDAs (Fig 2d and Fig 2e)
#1 and 2
plot(junco.prda, type="n", scaling=3, choices=c(1,2))
points(junco.prda, display="sites", pch=21, cex=2, col="gray32", scaling=3, bg=myCol) # the juncos
text(junco.prda, scaling=3, display="bp", cex=1, choices=c(1,2)) # the predictors

#3 and 4
plot(junco.prda, type="n", scaling=3, choices=c(3,4))
points(junco.prda, display="sites", pch=21, cex=2, col="gray32", scaling=3, choices=c(3,4), bg=myCol) # the juncos
text(junco.prda, scaling=3, display="bp", cex=1, choices=c(3,4)) # the predictors
