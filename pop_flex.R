#This script runs analyses to quantify differences among acclimated Junco populations as shown in: Stager M, Senner NR, Swanson DL, Carling MD, Grieves TJ, and Cheviron ZA. 2021. Temperature heterogeneity correlates with intraspecific variation in physiological flexibility in a small endotherm. Nature Communications, accepted

data<-read.csv("Source_Data_Fig3") #this dataset is included in the Source Data file for Figure 3 accompanying the paper

################################################
#Mantel Test

#sampling locations, in following order: J.h.shufeldti, J.h.aikeni, J.p.palliatus, J.h.dorsalis, J.h.thurberi
lats <- c(43.63, 44.26, 31.89, 34.43, 32.72)
lons <- c(-124.03, -103.81, -109.28, -111.40, -117.16) 


#pairwise geodesic (ellipsoid distances) between sampling points in km
library(geosphere)
geo_dist = matrix(NA, 5,5)
for (i in 1:5){
	for (j in 1:5) {
		geo_dist[i,j] = distGeo(c(lons[i],lats[i]),c(lons[j],lats[j]))/1000 #longitude, latitude
	}
}

#pairwise genetic distance: Fst/(1-Fst) using Weir's weighted theta, values reported in TableS4
#listed in order to match lats/longs above: J.h.shufeldti, J.h.aikeni, J.p.palliatus, J.h.dorsalis, J.h.thurberi
fst = matrix (c(c(NA,0.023624,0.051234,0.018739,0.018673), 
				c(0.023624,NA,0.051391,0.026072,0.026478), 
				c(0.051234,0.051391,NA,0.028971,0.041593), 
				c(0.018739,0.026072,0.028971,NA,0.01858), 
				c(0.018673,0.026478,0.041593,0.01858,NA)), nrow=5, ncol=5)

gen_dist = fst / (1-fst)

#Environmental distance 
#first extract climate variables from WorldClim dataset
library(raster)
clim=getData('worldclim', var='bio', res=2.5, download=TRUE)
r <- clim[[c(1:19)]]
values <- extract(r,SpatialPoints(data.frame(x=lons,y=lats), proj4string = r@crs))
#remove highly correlated variables at r > 0.7
library(psych)
pairs.panels(values, scale=T, cex=2) #bio1,4,5,6,7,9,11:15,17:19 redundant
pairs.panels(values[,c(2:3,8,10,16)], scale=T, cex=2) #no pairs show r > 0.7
env_dist = as.matrix(dist(values[,c(2:3,8,10,16)], method="euclidean", upper=TRUE, diag=TRUE))

#Partial Mantel test
library(vegan)
mantel.partial(gen_dist, env_dist, geo_dist) #IBE

detach("package:psych", unload=TRUE) #psych library cannot be loaded to use rescale function from arm

#Plot sampling sites (Fig 3a)
library(wesanderson)

ext = extent(-125,-95,30,50)
c2 = crop(data,ext)
c2.df = as.data.frame(c2, xy=TRUE)
States <- map_data("state")

ggplot() + 
	geom_raster(data = c2.df, aes(x = x, y = y, fill = bio7)) +
	coord_quickmap() + 
	xlab("lat") + ylab("long") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
	theme(panel.background=element_rect(fill="white")) + 
	theme(axis.line = element_line(colour = "black")) +
	geom_point(aes(y=lats, x=lons, size =2)) +
	geom_polygon( data=States, aes(x=long, y=lat, group=group),color="black", fill=NA) +
	scale_fill_gradientn(colours=wes_palette("Zissou1", 5,type="continuous"), na.value="white") + coord_cartesian(xlim=c(-125, -100), ylim = c(30, 50))

####################################################
#Add temp range data from WorldClim in deg C to dataset
data$range_temp[data$Population=="J. h. shufeldti"] <- values[1,7]/10
data$range_temp[data$Population=="J. h. aikeni"] <- values[2,7]/10
data$range_temp[data$Population=="J. p. palliatus"] <- values[3,7]/10
data$range_temp[data$Population=="J. h. dorsalis"] <- values[4,7]/10
data$range_temp[data$Population=="J. h. thurberi"] <- values[5,7]/10

#######################################################################
#influence of start time
data$Start_Time_PreAcc_Msum_Trial<-as.integer(format(as.POSIXct(data$Start_Time_PreAcc_Msum_Trial, format="%H:%M"), "%H"))
data$Start_Time_PostAcc_Msum_Trial<-as.integer(format(as.POSIXct(data$Start_Time_PostAcc_Msum_Trial, format="%H:%M"), "%H"))
t.test(Start_Time_PostAcc_Msum_Trial~Treatment, data) #different
summary(lm(PreAcc_Msum~Start_Time_PreAcc_Msum_Trial, data)) #no effect of start time on VO2
summary(lm(PreAcc_Msum~Start_Time_PostAcc_Msum_Trial, data)) #no effect of start time on VO2

#####################################
#Calculate difference in Msum with acclimation (post - pre-acclimation values)
data$delta_Msum = data$PostAcc_Msum - data$PreAcc_Msum

#coefficients of variation
trange <- c(data$range_temp[data$Population=="J. h. thurberi"][1], data$range_temp[data$Population=="J. h. shufeldti"][1], data$range_temp[data$Population=="J. p. palliatus"][1], data$range_temp[data$Population=="J. h. dorsalis"][1], data$range_temp[data$Population=="J. h. aikeni"][1])

#Cold birds CVs
cvs <- c(cv(data$delta_Msum[data$Population=="J. h. thurberi" & data$Treatment=="Cold"]), cv(data$delta_Msum[data$Population=="J. h. shufeldti" & data$Treatment=="Cold"]), cv(data$delta_Msum[data$Population=="J. p. palliatus" & data$Treatment=="Cold"]), cv(data$delta_Msum[data$Population=="J. h. dorsalis" & data$Treatment=="Cold"]), cv(data$delta_Msum[data$Population=="J. h. aikeni" & data$Treatment=="Cold"]))

summary(lm(cvs~trange))

#Control birds
cvs <- c(cv(data$delta_Msum[data$Population=="J. h. thurberi" & data$Treatment=="Control"]), cv(data$delta_Msum[data$Population=="J. h. shufeldti" & data$Treatment=="Control"]), cv(data$delta_Msum[data$Population=="J. p. palliatus" & data$Treatment=="Control"]), cv(data$delta_Msum[data$Population=="J. h. dorsalis" & data$Treatment=="Control"]), cv(data$delta_Msum[data$Population=="J. h. aikeni" & data$Treatment=="Control"]))

summary(lm(cvs~trange))
################################################
#Test for effect of temperature range on acclimation ability while controlling for population differentiation with MCMC GLMM

#First, use Fst values (Shown in Figure 3 or output from vcftools in accompanying script) to manually construct pairwise matrix for each individual
data_fst = matrix(nrow=nrow(data), ncol=nrow(data))
for (i in 1:nrow(data_fst)){
	for (j in 1:nrow(data_fst)){
		if (data$Population[i]==data$Population[j]) {data_fst[i,j] = 0}
		else if (data$Population[i]=="J. h. dorsalis" & data$Population[j]=="J. h. thurberi" | data$Population[i]=="J. h. thurberi" & data$Population[j]=="J. h. dorsalis") {data_fst[i,j] = 0.01858}
		else if (data$Population[i]=="J. h. dorsalis" & data$Population[j]=="J. h. shufeldti" | data$Population[i]=="J. h. shufeldti" & data$Population[j]=="J. h. dorsalis") {data_fst[i,j] = 0.018739}
		else if (data$Population[i]=="J. h. dorsalis" & data$Population[j]=="J. h. aikeni" | data$Population[i]=="J. h. aikeni" & data$Population[j]=="J. h. dorsalis") {data_fst[i,j] = 0.026072}
		else if (data$Population[i]=="J. h. dorsalis" & data$Population[j]=="J. p. palliatus" | data$Population[i]=="J. p. palliatus" & data$Population[j]=="J. h. dorsalis") {data_fst[i,j] = 0.028971}
		else if (data$Population[i]=="J. h. thurberi" & data$Population[j]=="J. h. shufeldti" | data$Population[i]=="J. h. shufeldti" & data$Population[j]=="J. h. thurberi") {data_fst[i,j] = 0.018673}
		else if (data$Population[i]=="J. h. thurberi" & data$Population[j]=="J. h. aikeni" | data$Population[i]=="J. h. aikeni" & data$Population[j]=="J. h. thurberi") {data_fst[i,j] = 0.026478}
		else if (data$Population[i]=="J. h. thurberi" & data$Population[j]=="J. p. palliatus" | data$Population[i]=="J. p. palliatus" & data$Population[j]=="J. h. thurberi") {data_fst[i,j] = 0.041593}
		else if (data$Population[i]=="J. h. shufeldti" & data$Population[j]=="J. h. aikeni" | data$Population[i]=="J. h. aikeni" & data$Population[j]=="J. h. shufeldti") {data_fst[i,j] = 0.023624}
		else if (data$Population[i]=="J. h. shufeldti" & data$Population[j]=="J. p. palliatus" | data$Population[i]=="J. p. palliatus" & data$Population[j]=="J. h. shufeldti") {data_fst[i,j] = 0.051234}
		else if (data$Population[i]=="J. h. aikeni" & data$Population[j]=="J. p. palliatus" | data$Population[i]=="J. p. palliatus" & data$Population[j]=="J. h. aikeni") {data_fst[i,j] = 0.051391}
	}
}

data$fst = data_fst #add this matrix to the dataframe

#Run MCMC GLMM model
library(MCMCglmm)
library(arm) #includes function rescale() for standardizing according to Gelman (2008)

#pre-treatment differences in mass?
summary(MCMCglmm(fixed=PreAcc_Mass_.g.~Treatment, random = ~fst, data = data, nitt=1000000, burnin=10000, thin=100, verbose=FALSE))

#pre-treatment differences in Msum?
summary(MCMCglmm(fixed=PreAcc_Msum~Treatment, random = ~fst, data = data, nitt=1000000, burnin=10000, thin=100, verbose=FALSE))

#pre-treatment differences in Mass by temp range?
summary(MCMCglmm(fixed=PreAcc_Mass_.g.~rescale(range_temp), random = ~fst, data = data, nitt=1000000, burnin=10000, thin=100, verbose=FALSE))

#pre-treatment differences in Msum by temp range?
summary(MCMCglmm(fixed=PreAcc_Msum~rescale(range_temp), random = ~fst, data = data, nitt=1000000, burnin=10000, thin=100, verbose=FALSE))

#effect of capture season on pre-acclimation Msum?
summary(MCMCglmm(fixed=PreAcc_Msum~rescale(PreAcc_Mass_.g.)+Capture.Season, random = ~fst, data = data, nitt=1000000, burnin=10000, thin=100, verbose=FALSE))

#VO2max at different temps in CTRL birds
t.test(pop$SMR1_VO2max[pop$Treatment=="Control"],pop$SMR2_VO2max[pop$Treatment=="Control"], paired=TRUE)

#Msum flexibility
mod1 <- (MCMCglmm(fixed=delta_Msum~rescale(PreAcc_Mass_.g.)+rescale(range_temp)*Treatment, random = ~fst, data = data, nitt=1000000, burnin=10000, thin=100, verbose=FALSE)) 
summary(mod1)
#plot(mod1)

#effect of capture season on flexibility in Msum?
summary(MCMCglmm(fixed=delta_Msum~rescale(PreAcc_Mass_.g.)+Capture.Season, random = ~fst, data = data, nitt=1000000, burnin=10000, thin=100, verbose=FALSE))

#Order populations from least to greatest temp range for plotting
data$Population = factor(data$Population, levels=c("J. h. thurberi","J. h. shufeldti","J. p. palliatus","J. h. dorsalis","J. h. aikeni"))
data$Treatment=factor(data$Treatment, levels=c("Control", "Cold"))

#Plot Fig3b
ggplot(data, aes(x=factor(Population), y=(delta_Msum), colour=Treatment))+ geom_hline(yintercept=0, color="grey") + geom_boxplot(outlier.shape = NA)+ geom_point(position=position_jitterdodge(dodge.width=.75, jitter.width=.2),alpha=0.5)  + xlab("Population") + ylab("Delta Thermogenic Capacity (ml O2/min)") + labs(fill="Legend")+ theme(panel.background=element_rect(fill="white"), plot.margin=unit(c(1,1,1,1),"cm"))+ theme(axis.title.y=element_text(margin=margin(0,20,0,0)), axis.title.x=element_text(margin=margin(20,0,0,0)))+theme(legend.position="none")+scale_colour_manual(values=c("indianred2", "dodgerblue3"))
