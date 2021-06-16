#This script runs analyses to quantify differences among Juncos caught and measured in situ, as shown in: Stager M, Senner NR, Swanson DL, Carling MD, Grieves TJ, and Cheviron ZA. 2021. Temperature heterogeneity correlates with intraspecific variation in physiological flexibility in a small endotherm. Nature Communications, accepted.

library(daymetr)
library(arm)
library(MuMIn)
library(wesanderson)

data = read.csv("Source_Data_Fig1.csv") #data can be found in the Source Data for Figure 1 accompanying the paper

#get elevation for each site
library(googleway)
data$elev <- google_elevation(data.frame(lat=data$Latitude, lon=data$Longitude),  key = "API_key_here")$results$elevation

#formatting dates
data$Capture.Date2 = as.Date(data$Capture.Date, "%m/%d/%y")
data$Capture.Date = format(data$Capture.Date2, "%d/%m/%Y")

#how long were individuals held?
data$captivity <- as.Date(data$Msum.Date, format="%m/%d/%y")-data$Capture.Date2
mean(data$captivity)

df = data.frame(unique(data$Latitude),unique(data$Longitude))
names(df) = c("Latitude", "Longitude")
for (i in 1:nrow(df)) {
	df$minDate[i] = as.character(min(data$Capture.Date2[data$Latitude==df$Latitude[i]], na.rm=TRUE))
	df$maxDate[i] = as.character(max(data$Capture.Date2[data$Latitude==df$Latitude[i]], na.rm=TRUE))
}
df$site = c(1:18) 

#Get interpolated weather data from Daymet for the year of and the year preceding measure
library(daymetr)
rm(tmin)
rm(tmax)
rm(prcp)
rm(dayl)
rm(wvp)
rm(srad)
for (i in 1:nrow(df)){
	clim <- download_daymet(site = df$site[i],
		lat = df$Latitude[i], 
		lon = df$Longitude[i], 
		start = as.numeric(substring(df$minDate[i], 1, 4))-1, 
		end = as.numeric(substring(df$maxDate[i], 1, 4)), 
		internal = TRUE, 
		simplify = TRUE)
	if (!exists("tmin")){tmin <- clim[grep("tmin", clim$measurement),]}
		else {tmin = rbind(tmin, clim[grep("tmin", clim$measurement),])}
	if (!exists("tmax")){tmax <- clim[grep("tmax", clim$measurement),]}
		else {tmax = rbind(tmax, clim[grep("tmax", clim$measurement),])}
	if (!exists("prcp")){prcp <- clim[grep("prcp", clim$measurement),]}
		else {prcp = rbind(prcp, clim[grep("prcp", clim$measurement),])}
	if (!exists("dayl")){dayl <- clim[grep("dayl", clim$measurement),]}
		else {dayl = rbind(dayl, clim[grep("dayl", clim$measurement),])}
	if (!exists("wvp")){wvp <- clim[grep("vp", clim$measurement),]}
		else {wvp = rbind(wvp, clim[grep("vp", clim$measurement),])}
	if (!exists("srad")){srad <- clim[grep("srad", clim$measurement),]}
		else {srad = rbind(srad, clim[grep("srad", clim$measurement),])}
}

tmin$Date = format(as.Date(as.numeric(tmin$yday)-1, origin=paste(tmin$year,"-01-01", sep="")), "%d/%m/%Y")
tmin$value = as.numeric(tmin$value)

tmax$Date = format(as.Date(as.numeric(tmax$yday)-1, origin=paste(tmax$year,"-01-01", sep="")), "%d/%m/%Y")
tmax$value = as.numeric(tmax$value)

prcp$Date = format(as.Date(as.numeric(prcp$yday)-1, origin=paste(prcp$year,"-01-01", sep="")), "%d/%m/%Y")
prcp$value = as.numeric(prcp$value)

dayl$Date = format(as.Date(as.numeric(dayl$yday)-1, origin=paste(dayl$year,"-01-01", sep="")), "%d/%m/%Y")
dayl$value = as.numeric(dayl$value)

wvp$Date = format(as.Date(as.numeric(wvp$yday)-1, origin=paste(wvp$year,"-01-01", sep="")), "%d/%m/%Y")
wvp$value = as.numeric(wvp$value)

srad$Date = format(as.Date(as.numeric(srad$yday)-1, origin=paste(srad$year,"-01-01", sep="")), "%d/%m/%Y")
srad$value = as.numeric(srad$value)

data2 <- data
data2$Taxon = as.factor(data2$Taxon)
data2$Taxon = relevel(data2$Taxon, ref="J. h. oreganus group")

#Add mean weather variables for 14 days preceding measurement
windows = data.frame(matrix(ncol=8,nrow=8)) #store AIC values in table (shown as Table S2)
names(windows)=c("window","tmin", "tmax", "prcp", "dayl", "wvp", "srad","range")

library(arm)
m = 1
for (k in 7:14) {
	for (i in 1:nrow(data2)){
			data2$tmin[i] <- mean(tmin$value[(which(tmin$Date == data2$Capture.Date[i] & tmin$latitude==data2$Latitude[i])-k):(which(tmin$Date == data2$Capture.Date[i] & tmin$latitude==data2$Latitude[i])-1)])
			data2$tmax[i] <- mean(tmax$value[(which(tmax$Date == data2$Capture.Date[i] & tmin$latitude==data2$Latitude[i])-k):(which(tmax$Date == data2$Capture.Date[i] & tmin$latitude==data2$Latitude[i])-1)])
			data2$range[i] <- mean(tmax$value[(which(tmax$Date == data2$Capture.Date[i] & tmax$latitude==data2$Latitude[i])-k):(which(tmax$Date == data2$Capture.Date[i] & tmax$latitude==data2$Latitude[i])-1)]-tmin$value[(which(tmin$Date == data2$Capture.Date[i] & tmin$latitude==data2$Latitude[i])-k):(which(tmin$Date == data2$Capture.Date[i] & tmin$latitude==data2$Latitude[i])-1)])
			data2$prcp[i] <- mean(prcp$value[(which(prcp$Date == data2$Capture.Date[i] & tmin$latitude==data2$Latitude[i])-k):(which(prcp$Date == data2$Capture.Date[i] & tmin$latitude==data2$Latitude[i])-1)])
			data2$dayl[i] <- mean(dayl$value[(which(dayl$Date == data2$Capture.Date[i] & tmin$latitude==data2$Latitude[i])-k):(which(dayl$Date == data2$Capture.Date[i] & tmin$latitude==data2$Latitude[i])-1)])
			data2$wvp[i] <- mean(wvp$value[(which(wvp$Date == data2$Capture.Date[i] & tmin$latitude==data2$Latitude[i])-k):(which(wvp$Date == data2$Capture.Date[i] & tmin$latitude==data2$Latitude[i])-1)])
			data2$srad[i] <- mean(srad$value[(which(srad$Date == data2$Capture.Date[i] & tmin$latitude==data2$Latitude[i])-k):(which(srad$Date == data2$Capture.Date[i] & tmin$latitude==data2$Latitude[i])-1)])
	}
	windows$window[m]<-paste(1,':', k)
	windows$tmin[m]<-AIC(lm(Msum~rescale(Mass) + Taxon * rescale(tmin), data2)) 
	windows$tmax[m]<-AIC(lm(Msum~rescale(Mass) + Taxon * rescale(tmax), data2)) 
	windows$prcp[m]<-AIC(lm(Msum~rescale(Mass) + Taxon * rescale(prcp), data2)) 
	windows$dayl[m]<-AIC(lm(Msum~rescale(Mass) + Taxon * rescale(dayl), data2)) 
	windows$wvp[m]<-AIC(lm(Msum~rescale(Mass) + Taxon * rescale(wvp), data2)) 
	windows$range[m]<-AIC(lm(Msum~rescale(Mass) + Taxon * rescale(range), data2))
	windows$srad[m]<-AIC(lm(Msum~rescale(Mass) + Taxon * rescale(srad), data2))
	m = m+1 
} 

AIC(lm(Msum~rescale(Mass) + Taxon, data2)) #null
AIC(lm(Msum~rescale(Mass) + Taxon + rescale(elev), data2)) 	

#best window from above
j = 1
k = 8
for (i in 1:nrow(data2)){
		data2$range[i] <- mean(tmax$value[(which(tmax$Date == data2$Capture.Date[i] & tmax$latitude==data2$Latitude[i])-k):(which(tmax$Date == data2$Capture.Date[i] & tmax$latitude==data2$Latitude[i])-j)]-tmin$value[(which(tmin$Date == data2$Capture.Date[i] & tmin$latitude==data2$Latitude[i])-k):(which(tmin$Date == data2$Capture.Date[i] & tmin$latitude==data2$Latitude[i])-j)])
	}

cor(data2$tmin, data2$tmax) #0.96

data2$Taxon = relevel(data2$Taxon, ref="J. h. oreganus group")
mod1 <- lm(Msum~rescale(Mass) + Taxon * rescale(range), data2)#**** 
summary(mod1)

#significance of interaction term
AIC(mod1)-AIC(lm(Msum~rescale(Mass) + Taxon + rescale(range), data2))

#significance of Season 
AIC(mod1)-AIC(lm(Msum~rescale(Mass) + Taxon * rescale(range) + Season, data2))



#plotting model results
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

df_caniceps = data.frame(range = seq(from = min(data2$range), to = max(data2$range), by = .25), Mass = mean(data2$Mass[data2$Taxon=="J. h. caniceps"]), Taxon = "J. h. caniceps")
df_mearnsi = data.frame(range = seq(from = min(data2$range), to = max(data2$range), by = .25), Mass = mean(data2$Mass[data2$Taxon=="J. h. mearnsi"]), Taxon = "J. h. mearnsi")
df_oreganus= data.frame(range = seq(from = min(data2$range), to = max(data2$range), by = .25), Mass = mean(data2$Mass[data2$Taxon=="J. h. oreganus group"]), Taxon = "J. h. oreganus group")
df_hyemalis= data.frame(range = seq(from = min(data2$range), to = max(data2$range), by = .25), Mass = mean(data2$Mass[data2$Taxon=="J. h. hyemalis"]), Taxon = "J. h. hyemalis")
df_palliatus = data.frame(range = seq(from = min(data2$range), to = max(data2$range), by = .25), Mass = mean(data2$Mass[data2$Taxon=="J. p. palliatus"]), Taxon = "J. p. palliatus")

df_caniceps$prediction<-predict(lm(Msum~Mass + Taxon * range, data2), newdata=df_caniceps, type="response")
df_mearnsi$prediction<-predict(lm(Msum~Mass + Taxon * range, data2), newdata=df_mearnsi, type="response")
df_oreganus$prediction<-predict(lm(Msum~Mass + Taxon * range, data2), newdata=df_oreganus, type="response")
df_hyemalis$prediction<-predict(lm(Msum~Mass + Taxon * range, data2), newdata=df_hyemalis, type="response")
df_palliatus$prediction<-predict(lm(Msum~Mass + Taxon * range, data2), newdata=df_palliatus, type="response")


#Plot Figure 1
library(wesanderson); library(ggplot2)
pal <- wes_palette(5, name = "Darjeeling1", type = "continuous")
data2$Taxon = factor(data2$Taxon, levels=c("J. p. palliatus","J. h. caniceps","J. h. mearnsi","J. h. oreganus group","J. h. hyemalis"))
df$pal = c(3,3,3,3,5,5,5,5,1,4,4,2,1,1,2,2,2,2) #add colors to df


#Fig. 1a - Map
library(mapdata); library(maptools)
usa <- map_data("state")
canada <- map_data("worldHires", "Canada")
mexico <- map_data("worldHires", "Mexico")

sites <- ggplot() + 
	geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = "white", color="black") + 	
	geom_polygon(data = canada, aes(x=long, y = lat, group = group), fill = "white", color="black") + 
	geom_polygon(data = mexico, aes(x=long, y = lat, group = group), fill = "white", color="black") + 
	coord_fixed(xlim = c(-123, -69),  ylim = c(25, 49), ratio = 1.7) + 
	geom_point(data=df[c(1:10,12:15),], aes(y=Latitude, x=Longitude), shape=21, fill=pal[df$pal[c(1:10,12:15)]],colour="black", size=7) +
	theme(plot.margin=unit(c(0,0,0,1),"cm")) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
	theme(axis.line = element_line(colour = "black")) +
	theme(axis.title.y=element_text(margin=margin(0,10,0,0)), axis.title.x=element_text(margin=margin(10,0,0,0)))

#Fig 1b - J. h. caniceps data
caniceps <- ggplot() + 
	xlim(7,23) + ylim(1,11) +
	xlab(expression(paste("Daily Temperature Range "( degree*C)))) + ylab("Msum") +
	geom_point(data = data2[data2$Taxon=="J. h. caniceps" & data2$Season=="Breeding",], aes(y=Msum, x=range), color=pal[5], shape=19, size = 2.5) +
	geom_point(data = data2[data2$Taxon=="J. h. caniceps" & data2$Season=="Non-breeding",], aes(y=Msum, x=range), shape=21, fill="white", color=pal[5], size = 2.5) + 
	geom_abline(intercept = coef(lm(prediction~range, df_caniceps))[1], slope = coef(lm(prediction~range, df_caniceps))[2],  colour=pal[5], linetype=2, size = 1) + 		
	theme(panel.background=element_rect(fill="white"), plot.margin=unit(c(1,1,1,1),"cm")) + 
	theme(axis.line = element_line(colour = "black")) +
	theme(axis.title.y=element_text(margin=margin(0,10,0,0)), axis.title.x=element_text(margin=margin(10,0,0,0)))
	
#Fig 1c - J. h. hymealis data
hyemalis <- ggplot() + xlim(7,23) + ylim(1,11) +
	xlab(expression(paste("Daily Temperature Range "( degree*C)))) + ylab("Msum") +
	geom_point(data = data2[data2$Taxon=="J. h. hyemalis" & data2$Season=="Breeding",], aes(y=Msum, x=range), color=pal[1], shape=19, size = 2.5) +
	geom_point(data = data2[data2$Taxon=="J. h. hyemalis" & data2$Season=="Non-breeding",], aes(y=Msum, x=range), shape=21, fill="white", color=pal[1], size = 2.5) + 
	geom_abline(intercept = coef(lm(prediction~range, df_hyemalis))[1], slope = coef(lm(prediction~range, df_hyemalis))[2],  colour=pal[1], linetype=2, size = 1) + 		
	theme(panel.background=element_rect(fill="white"), plot.margin=unit(c(1,1,1,1),"cm")) + 
	theme(axis.line = element_line(colour = "black")) +
	theme(axis.title.y=element_text(margin=margin(0,10,0,0)), axis.title.x=element_text(margin=margin(10,0,0,0)))+
	theme(legend.position="none")
	
#Fig 1d - J. h. mearnsi data 
mearnsi <- ggplot() + xlim(7,23) + ylim(1,11) +
	xlab(expression(paste("Daily Temperature Range "( degree*C)))) + ylab("Msum") +
	geom_point(data = data2[data2$Taxon=="J. h. mearnsi" & data2$Season=="Breeding",], aes(y=Msum, x=range), color=pal[2], shape=19, size = 2.5) +
	geom_point(data = data2[data2$Taxon=="J. h. mearnsi" & data2$Season=="Non-breeding",], aes(y=Msum, x=range), shape=21, fill="white", color=pal[2], size = 2.5) + 
	geom_abline(intercept = coef(lm(prediction~range, df_mearnsi))[1], slope = coef(lm(prediction~range, df_mearnsi))[2],  colour=pal[2], linetype=2, size = 1) + 		
	theme(panel.background=element_rect(fill="white"), plot.margin=unit(c(1,1,1,1),"cm")) + 
	theme(axis.line = element_line(colour = "black")) +
	theme(axis.title.y=element_text(margin=margin(0,10,0,0)), axis.title.x=element_text(margin=margin(10,0,0,0)))+
	theme(legend.position="none")

#Fig 1e - J. h. oreganus group data
oreganus <- ggplot() + xlim(7,23) + ylim(1,11) +
	xlab(expression(paste("Daily Temperature Range "( degree*C)))) + ylab("Msum") +
	geom_point(data = data2[data2$Taxon=="J. h. oreganus group" & data2$Season=="Breeding",], aes(y=Msum, x=range), color=pal[4], shape=19, size = 2.5) +
	geom_point(data = data2[data2$Taxon=="J. h. oreganus group" & data2$Season=="Non-breeding",], aes(y=Msum, x=range), shape=21, fill="white", color=pal[4], size = 2.5) + 
	geom_abline(intercept = coef(lm(prediction~range, df_oreganus))[1], slope = coef(lm(prediction~range, df_oreganus))[2],  colour=pal[4], linetype=2, size = 1) + 		
	theme(panel.background=element_rect(fill="white"), plot.margin=unit(c(1,1,1,1),"cm")) + 
	theme(axis.line = element_line(colour = "black")) +
	theme(axis.title.y=element_text(margin=margin(0,10,0,0)), axis.title.x=element_text(margin=margin(10,0,0,0)))+
	theme(legend.position="none")

#Fig 1f - J. p. palliatus data
palliatus <- ggplot() + xlim(7,23) + ylim(1,11) +
	xlab(expression(paste("Daily Temperature Range "( degree*C)))) + ylab("Msum") +
	geom_point(data = data2[data2$Taxon=="J. p. palliatus" & data2$Season=="Breeding",], aes(y=Msum, x=range), color=pal[3], shape=19, size = 2.5) +
	geom_point(data = data2[data2$Taxon=="J. p. palliatus" & data2$Season=="Non-breeding",], aes(y=Msum, x=range), shape=21, fill="white", color=pal[3], size = 2.5) + 
	geom_abline(intercept = coef(lm(prediction~range, df_palliatus))[1], slope = coef(lm(prediction~range, df_palliatus))[2],  colour=pal[3], linetype=2, size = 1) + 		
	theme(panel.background=element_rect(fill="white"), plot.margin=unit(c(1,1,1,1),"cm")) + 
	theme(axis.line = element_line(colour = "black")) +
	theme(axis.title.y=element_text(margin=margin(0,10,0,0)), axis.title.x=element_text(margin=margin(10,0,0,0)))+
	theme(legend.position="none")

multiplot(sites,mearnsi,caniceps,oreganus,hyemalis,palliatus, cols=3) #plot six-paneled figure


#Reaction norms vs. temperature range
palliatus_slope <- (df_palliatus$prediction[2]-df_palliatus$prediction[1])/.25
caniceps_slope <- (df_caniceps$prediction[2]-df_caniceps$prediction[1])/.25
mearnsi_slope <- (df_mearnsi$prediction[2]-df_mearnsi$prediction[1])/.25
oreganus_slope <- (df_oreganus$prediction[2]-df_oreganus$prediction[1])/.25
hyemalis_slope <- (df_hyemalis$prediction[2]-df_hyemalis$prediction[1])/.25

slopes <-c(palliatus_slope, caniceps_slope, mearnsi_slope, oreganus_slope, hyemalis_slope)

library(raster)
data=getData('worldclim', var='bio', res=2.5, download=TRUE) #bioclim data
r <- data[[7]] #Annual Temp Range
names(r) <- "Annual.Range"
coords <- data.frame(long=data2$Longitude,lat=data2$Latitude)
points <- SpatialPoints(coords, proj4string = r@crs)
values <- extract(r,points) #extract data for sample locales
clim <- cbind.data.frame(coordinates(points),values)
data2$Annual.Range <- clim$values/10 #degrees C

Annual.Range <- c(mean(data2$Annual.Range[data2$Taxon=="J. p. palliatus"]), mean(data2$Annual.Range[data2$Taxon=="J. h. caniceps"]), mean(data2$Annual.Range[data2$Taxon=="J. h. mearnsi"]), mean(data2$Annual.Range[data2$Taxon=="J. h. oreganus group"]), mean(data2$Annual.Range[data2$Taxon=="J. h. hyemalis"]))

summary(lm(slopes~Annual.Range))
summary(lm(slopes[-3]~Annual.Range[-3])) #w/o J. h. mearnsi

#Plot Figure S1
ggplot() + 
	xlab(expression(paste("Annual Temperature Range "( degree*C)))) + ylab("Flexibility in Msum (slope)") +
	geom_abline(intercept = coef(lm(slopes[-3]~Annual.Range[-3]))[1], slope = coef(lm(slopes[-3]~Annual.Range[-3]))[2],   linetype=2, size = 1, col="gray") + 
geom_point(aes(x=Annual.Range, y=slopes), color=pal[c(3,5,2,4,1)], shape=19, size = 4) +
	theme(panel.background=element_rect(fill="white"), plot.margin=unit(c(1,1,1,1),"cm")) + 
	theme(axis.line = element_line(colour = "black")) +
	theme(axis.title.y=element_text(margin=margin(0,10,0,0)), axis.title.x=element_text(margin=margin(10,0,0,0)))+
	theme(legend.position="none")
