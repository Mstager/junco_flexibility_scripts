#This script runs a redundancy analysis to identify genotype-environment associations, as shown in: Stager M, Senner NR, Swanson DL, Carling MD, Grieves TJ, and Cheviron ZA. 2020. Temperature heterogeneity correlates with intraspecific variation in physiological flexibility in a small endotherm. bioRxiv.

#This workflow has been modified from Forester et al. 2018's tutorial (available at https://popgen.nescent.org/2018-03-27_RDA_GEA.html) to incorporate variance partitioning

setwd("/Users/Maria/Google_Drive/Manuscripts/Dissertation/Ch4/Submission/")
#################################################
library(adegenet)
library(pegas)
library(vegan)
library(factoextra)

SNPS <- read.PLINK("github/Junco_SNPs.raw") 
gen <- as.data.frame(SNPS)
dim(gen)

geo<-read.csv("Supplement/TableS2.csv") #georeferenced samples
geo <- geo[geo$Catalog_number!="6056" & geo$Catalog_number!="110508" & geo$Catalog_number!="117690" & geo$Catalog_number!="118252" & geo$Catalog_number!="229176" & geo$Catalog_number!="231996" & geo$Catalog_number!="648103" & geo$Catalog_number!="648671" & geo$Catalog_number!="232004" & geo$Catalog_number!="183387" & geo$Catalog_number!="232018",] #remove 11 individuals that were filtered during RAD processing

sum(is.na(gen)) #total number of sites to impute
sum(is.na(gen)) / (dim(gen)[1]*dim(gen)[2]) #percentage of total


#impute using most common genotype for that population
gen.imp <- apply(gen[geo$Species=="hyemalis",], 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
gen.imp <- rbind(gen.imp, apply(gen[geo$Species=="insularis",], 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x)))))))
gen.imp <- rbind(gen.imp, apply(gen[geo$Species=="phaeonotus",], 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x)))))))
gen.imp <- rbind(gen.imp, apply(gen[geo$Species=="vulcani",], 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x)))))))

gen.imp = gen.imp[row.names(gen.imp)[match(row.names(gen), row.names(gen.imp))],]
identical(row.names(gen.imp[row.names(gen.imp)[match(row.names(gen), row.names(gen.imp))],]), geo$Catalog_number) #make sure files are in the same order

sum(is.na(gen.imp)) # No NAs

#PCA
pca_all <- dudi.pca(gen.imp,scannf=FALSE,scale=FALSE, nf=4)
get_eig(pca_all)[1:4,]

#Plotting
temp <- as.integer(as.factor(paste(geo$Species, geo$Subspecies)))
myCol <- c("green","coral2","cadetblue","orange","cyan","dodgerblue","deeppink","firebrick4","red","goldenrod1","darkseagreen1", "plum3","purple4","sienna4", "blue", "tan1", "darkgray", "yellowgreen","thistle3", "yellow2", "darkolivegreen")[temp] #alphabetical order

#plot PC axes 1 and 2
plot(pca_all$li[,1:2], col =myCol, pch=19)
legend("topright",col=unique(myCol),pch=19, legend=unique(paste(geo$Species, geo$Subspecies)), bty="n", cex=.7)

#plot PC axes 3 and 4
plot(pca_all$li[,c(3,4)], col =myCol, pch=19)
legend("topright",col=unique(myCol),pch=19, legend=unique(paste(geo$Species, geo$Subspecies)), bty="n", cex=.7)


#####################
#get environmental data
library(raster)

data=getData('worldclim', var='bio', res=2.5, download=TRUE) #bioclim data
r <- data[[c(1:19)]] #using all bioclim variables
names(r) <- c("Mean.Temp", "Diurnal.Range", "Isothermality", "Temp.Seasonality", "Max.Temp", "Min.Temp", "Annual.Range", "Temp.WettestQ", "Temp.DriestQ", "Temp.WarmestQ", "Temp.ColdestQ", "Annual.Precip", "Max.Precip", "Min.Precip", "Precip.Seasonality", "Precip.WettestQ", "Precip.DriestQ", "Precip.WarmestQ", "Precip.ColdestQ") #bioclim variable names
points <- SpatialPoints(data.frame(long=geo$Longitude,lat=geo$Latitude), proj4string = r@crs)
values <- extract(r,points) #extract data for sample locales
clim <- cbind.data.frame(coordinates(points),values)
clim $sampleID <- geo$Catalog_number #add catalog number
identical(rownames(gen.imp), clim$sampleID) #in same order as SNP data?

#scale variables
library(arm)
clim = data.frame(apply(clim, 2, rescale))

#removing environmental variables that are strongly correlated (r > 0.7)
library(psych)
pairs.panels(clim[,3:21], scale=T, cex=2.5)
pred <- subset(clim, select=-c(long,lat,sampleID,Mean.Temp, Temp.ColdestQ, Min.Temp, Annual.Precip, Temp.Seasonality, Max.Precip, Precip.DriestQ, Isothermality, Precip.WettestQ, Max.Temp, Precip.Seasonality)) #in order of removal, all r < 0.7
pairs.panels(pred, scale=T, cex=2.5)


#run the  rda
junco.rda <- rda(gen.imp ~ Diurnal.Range + Annual.Range + Temp.WettestQ + Temp.DriestQ + Temp.WarmestQ + Min.Precip + Precip.WarmestQ + Precip.ColdestQ, data = pred, scale=T)
junco.rda
RsquareAdj(junco.rda) #adjusted R-squared for number of predictors
summary(eigenvals(junco.rda, model = "constrained"))
vif.cca(junco.rda) 

#evaluate significance of model using permutation
signif.full.rda <- anova.cca(junco.rda, parallel=getOption("mc.cores")) # default is permutation=999, takes a while to run
signif.full.rda 

#evaluate significance of axes using permutation
signif.axis.rda <- anova.cca(junco.rda, by="axis", parallel=getOption("mc.cores"),  permutations=how(nperm=99)) #takes quite a while to run, runs out of memory at 999 permutations
signif.axis.rda 

#Variance partitioning
varpart(gen.imp, ~Diurnal.Range, ~ Annual.Range + Temp.WettestQ + Temp.DriestQ + Temp.WarmestQ + Min.Precip + Precip.WarmestQ + Precip.ColdestQ, data=pred) #amount of variance explained by Diurnal Range, represented by [a]

varpart(gen.imp, ~Annual.Range, ~ Diurnal.Range + Temp.WettestQ + Temp.DriestQ + Temp.WarmestQ + Min.Precip + Precip.WarmestQ + Precip.ColdestQ, data=pred)

varpart(gen.imp, ~ Temp.WettestQ, ~Diurnal.Range + Annual.Range + Temp.DriestQ + Temp.WarmestQ + Min.Precip + Precip.WarmestQ + Precip.ColdestQ, data=pred)

varpart(gen.imp, ~ Temp.DriestQ, ~Diurnal.Range + Annual.Range + Temp.WettestQ + Temp.WarmestQ + Min.Precip + Precip.WarmestQ + Precip.ColdestQ, data=pred)

varpart(gen.imp, ~ Temp.WarmestQ, ~Diurnal.Range + Annual.Range + Temp.WettestQ + Temp.DriestQ + Min.Precip + Precip.WarmestQ + Precip.ColdestQ, data=pred)

varpart(gen.imp, ~Min.Precip, ~Diurnal.Range + Annual.Range + Temp.WettestQ + Temp.DriestQ + Temp.WarmestQ + Precip.WarmestQ + Precip.ColdestQ, data=pred)

varpart(gen.imp, ~Precip.WarmestQ, ~Diurnal.Range + Annual.Range + Temp.WettestQ + Temp.DriestQ+  Temp.WarmestQ + Min.Precip + Precip.ColdestQ, data=pred)

varpart(gen.imp, ~Precip.ColdestQ, ~Diurnal.Range + Annual.Range + Temp.WettestQ + Temp.DriestQ+  Temp.WarmestQ + Min.Precip +Precip.WarmestQ, data=pred)


subsp<- as.factor(geo$Subspecies)

#Plot RDA axes 1 and 2
plot(junco.rda, type="n", scaling=3, xlim=c(-7,6))
points(junco.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=myCol) # the juncos
text(junco.rda, scaling=3, display="bp", col="gray32", cex=.8) # the predictors
#legend("topleft", legend=unique(paste(geo$Species, geo$Subspecies)), bty="n", col="gray32", pch=21, cex=.75, pt.bg=unique(myCol))

#Plot RDA axes 3 and 4
plot(junco.rda, type="n", scaling=3, choices=c(3,4), xlim=c(-8,4))
points(junco.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, choices=c(3,4), bg=myCol) # the juncos
text(junco.rda, scaling=3, display="bp", col="gray32", cex=.8, choices=c(3,4)) # the predictors
#legend("bottomleft", legend=unique(paste(geo$Species, geo$Subspecies)), bty="n", col="gray32", pch=21, cex=.75, pt.bg=unique(myCol))



####################################
#run the partial rda

#use PC axes 1 and 2 to condition RDA on background structure
pred$PC1 <- as.vector(pca_all$li[,1])
pred$PC2 <- pca_all$li[,2]

junco.prda <- rda(gen.imp ~ Diurnal.Range + Annual.Range + Temp.WettestQ + Temp.DriestQ + Temp.WarmestQ + Min.Precip + Precip.WarmestQ + Precip.ColdestQ + Condition(PC1+PC2), data=pred, scale=T) #same environmental variables as above
junco.prda
RsquareAdj(junco.prda) #adjusted R-squared for number of predictors
summary(eigenvals(junco.prda, model = "constrained"))
vif.cca(junco.prda) 

#evaluate significance of model using permutation
signif.full.prda <- anova.cca(junco.prda, parallel=getOption("mc.cores")) # default is permutation=999
signif.full.prda

#evaluate significance of axes using permutation
signif.axis.prda <- anova.cca(junco.prda, by="axis", parallel=getOption("mc.cores"),  permutations=how(nperm=99))
signif.axis.prda 

#variance partitioning
varpart(gen.imp, ~Diurnal.Range, ~ Annual.Range + Temp.WettestQ + Temp.DriestQ + Temp.WarmestQ + Min.Precip + Precip.WarmestQ + Precip.ColdestQ, ~PC1 + PC2, data=pred) #amount of variance explained by Diurnal Range, represented by [a]

varpart(gen.imp, ~Annual.Range, ~ Diurnal.Range + Temp.WettestQ + Temp.DriestQ + Temp.WarmestQ + Min.Precip + Precip.WarmestQ + Precip.ColdestQ, ~PC1 + PC2, data=pred)

varpart(gen.imp, ~ Temp.WettestQ, ~Diurnal.Range + Annual.Range + Temp.DriestQ + Temp.WarmestQ + Min.Precip + Precip.WarmestQ + Precip.ColdestQ, ~PC1 + PC2, data=pred)

varpart(gen.imp, ~ Temp.DriestQ, ~Diurnal.Range + Annual.Range + Temp.WettestQ + Temp.WarmestQ + Min.Precip + Precip.WarmestQ + Precip.ColdestQ, ~PC1 + PC2, data=pred)

varpart(gen.imp, ~ Temp.WarmestQ, ~Diurnal.Range + Annual.Range + Temp.WettestQ + Temp.DriestQ + Min.Precip + Precip.WarmestQ + Precip.ColdestQ, ~PC1 + PC2, data=pred)

varpart(gen.imp, ~Min.Precip, ~Diurnal.Range + Annual.Range + Temp.WettestQ + Temp.DriestQ + Temp.WarmestQ + Precip.WarmestQ + Precip.ColdestQ, ~PC1 + PC2, data=pred)

varpart(gen.imp, ~Precip.WarmestQ, ~Diurnal.Range + Annual.Range + Temp.WettestQ + Temp.DriestQ+  Temp.WarmestQ + Min.Precip + Precip.ColdestQ, ~PC1 + PC2, data=pred)

varpart(gen.imp, ~Precip.ColdestQ, ~Diurnal.Range + Annual.Range + Temp.WettestQ + Temp.DriestQ+  Temp.WarmestQ + Min.Precip  + Precip.WarmestQ, ~PC1 + PC2, data=pred)

varpart(gen.imp, ~PC1+PC2, ~Diurnal.Range + Annual.Range + Temp.WettestQ + Temp.DriestQ+  Temp.WarmestQ + Min.Precip  + Precip.WarmestQ + Precip.ColdestQ, data=pred)


#Plot RDA axes 1 and 2 from partial RDA
plot(junco.prda, type="n", scaling=3, xlim=c(-7,6))
points(junco.prda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=myCol) # the juncos
text(junco.prda, scaling=3, display="bp", col="gray32", cex=.8) # the predictors
#legend("topleft", legend=unique(paste(geo$Species, geo$Subspecies)), bty="n", col="gray32", pch=21, cex=.75, pt.bg=unique(myCol))

#Plot RDA axes 3 and 4 from partial RDA
plot(junco.prda, type="n", scaling=3, choices=c(2,3), xlim=c(-8,4))
points(junco.prda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, choices=c(2,3), bg=myCol) # the juncos
text(junco.rda, scaling=3, display="bp", col="gray32", cex=.8, choices=c(2,3)) # the predictors
#legend("bottomleft", legend=unique(paste(geo$Species, geo$Subspecies)), bty="n", col="gray32", pch=21, cex=.75, pt.bg=unique(myCol))

###################
#identify SNP candidates
load.prda <- scores(junco.prda, choices=c(1:3), display="species")  # Species scores for the first three constrained axes
par(mfrow=c(2,2))
hist(load.prda[,1], main="Loadings on RDA1")
hist(load.prda[,2], main="Loadings on RDA2")
hist(load.prda[,3], main="Loadings on RDA3") 


#identify SNPS in the tails of those distributions
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

cand1 <- outliers(load.prda[,1],3) 
cand2 <- outliers(load.prda[,2],3) 
cand3 <- outliers(load.prda[,3],3)
ncand <- length(cand1) + length(cand2) + length(cand3)
ncand  # Total number of outlier SNPs

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading") 

cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

#add correlations w/ env variables
foo <- matrix(nrow=(ncand), ncol=8)  # 8 columns for 8 predictors
colnames(foo) <- c("Diurnal.Range","Annual.Range","Temp.WettestQ","Temp.DriestQ", "Temp.WarmestQ","Min.Precip","Precip.WarmestQ", "Precip.ColdestQ")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- gen.imp[,nam]
  foo[i,] <- apply(pred[,1:8],2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)

length(cand$snp[duplicated(cand$snp)])  
foo <- cbind(cand$axis, duplicated(cand$snp)) 
table(foo[foo[,1]==1,2]) # no duplicates on axis 1
table(foo[foo[,1]==2,2]) #  0 duplicates on axis 2
table(foo[foo[,1]==3,2]) #  4 duplicates on axis 3
cand <- cand[!duplicated(cand$snp),] # remove duplicate detections
dim(cand) # SNPs remaining

#which predictor is each SNP most strongly correlated w/
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,12] <- names(which.max(abs(bar[4:11]))) # gives the variable
  cand[i,13] <- max(abs(bar[4:11]))              # gives the correlation
}

colnames(cand)[12] <- "predictor"
colnames(cand)[13] <- "correlation"

table(cand$predictor) 