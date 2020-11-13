#This script runs analyses to quantify differences among Juncos caught and measured in situ, as shown in: Stager M, Senner NR, Swanson DL, Carling MD, Grieves TJ, and Cheviron ZA. 2020. Temperature heterogeneity correlates with intraspecific variation in physiological flexibility in a small endotherm. bioRxiv.

library(daymetr)
library(arm)
library(MuMIn)
library(wesanderson)

data = read.csv("TableS1.csv")

#get elevation for each site
library(googleway)
data$elev <- google_elevation(data.frame(lat=data$Latitude, lon=data$Longitude),  key = "API_KEY_here")$results$elevation

#formatting dates
data$Capture.Date2 = as.Date(data$Capture.Date, "%m/%d/%y")
data$Capture.Date = format(data$Capture.Date2, "%d/%m/%Y")


df = data.frame(unique(data$Latitude),unique(data$Longitude))
names(df) = c("Latitude", "Longitude")
for (i in 1:nrow(df)) {
	df$minDate[i] = as.character(min(data$Capture.Date2[data$Latitude==df$Latitude[i]], na.rm=TRUE))
	df$maxDate[i] = as.character(max(data$Capture.Date2[data$Latitude==df$Latitude[i]], na.rm=TRUE))
}
df$site = c(1:18) 

#Get interpolated weather data from Daymet for the year of and the year preceding measure
rm(tmin)
rm(tmax)
rm(prcp)
rm(dayl)
rm(wvp)
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


data2 <- data
data2$climate <- 1
data2$Morphotype = as.factor(data2$Morphotype)
data2$Morphotype = relevel(data2$Morphotype, ref="Yellow-eyed")

#Add mean weather variables over acclimatization window (in days)

window = 8 #tested all windows from 7 to 14 d

for (i in 1:nrow(data2)){
	data2$tmin[i] <- mean(tmin$value[(which(tmin$Date == data2$Capture.Date[i] & tmin$latitude==data2$Latitude[i])-window):(which(tmin$Date == data2$Capture.Date[i] & tmin$latitude==data2$Latitude[i]))])
	data2$tmax[i] <- mean(tmax$value[(which(tmax$Date == data2$Capture.Date[i] & tmin$latitude==data2$Latitude[i])-window):(which(tmax$Date == data2$Capture.Date[i] & tmin$latitude==data2$Latitude[i]))])
	data2$range[i] <- mean(tmax$value[(which(tmax$Date == data2$Capture.Date[i] & tmax$latitude==data2$Latitude[i])-window):(which(tmax$Date == data2$Capture.Date[i] & tmax$latitude==data2$Latitude[i]))]-tmin$value[(which(tmin$Date == data2$Capture.Date[i] & tmin$latitude==data2$Latitude[i])-window):(which(tmin$Date == data2$Capture.Date[i] & tmin$latitude==data2$Latitude[i]))])
	data2$prcp[i] <- mean(prcp$value[(which(prcp$Date == data2$Capture.Date[i] & tmin$latitude==data2$Latitude[i])-window):(which(prcp$Date == data2$Capture.Date[i] & tmin$latitude==data2$Latitude[i]))])
	data2$dayl[i] <- mean(dayl$value[(which(dayl$Date == data2$Capture.Date[i] & tmin$latitude==data2$Latitude[i])-window):(which(dayl$Date == data2$Capture.Date[i] & tmin$latitude==data2$Latitude[i]))])
	data2$wvp[i] <- mean(wvp$value[(which(wvp$Date == data2$Capture.Date[i] & tmin$latitude==data2$Latitude[i])-window):(which(wvp$Date == data2$Capture.Date[i] & tmin$latitude==data2$Latitude[i]))])
	} 


AICc(lm(Msum_.ml.O2.min.~rescale(Mass_.g.) + Morphotype, data2)) #null
AICc(lm(Msum_.ml.O2.min.~rescale(Mass_.g.) + Morphotype + rescale(elev), data2)) 	
AICc(lm(Msum_.ml.O2.min.~rescale(Mass_.g.) + Morphotype + rescale(tmin), data2)) 
AICc(lm(Msum_.ml.O2.min.~rescale(Mass_.g.) + Morphotype + rescale(tmax), data2)) 
AICc(lm(Msum_.ml.O2.min.~rescale(Mass_.g.) + Morphotype + rescale(prcp), data2)) 
AICc(lm(Msum_.ml.O2.min.~rescale(Mass_.g.) + Morphotype + rescale(dayl), data2)) 
AICc(lm(Msum_.ml.O2.min.~rescale(Mass_.g.) + Morphotype + rescale(wvp), data2)) 
AICc(lm(Msum_.ml.O2.min.~rescale(Mass_.g.) + Morphotype + rescale(range), data2)) 

#Test interaction term for most well-supported model
summary(lm(Msum_.ml.O2.min.~rescale(Mass_.g.) + Morphotype * rescale(range), data2)) #**** 
data2$Morphotype = relevel(data2$Morphotype, ref="Yellow-eyed") #change reference population and rerun above


#plotting model results
df_GH = data.frame(range = seq(from = min(data2$range), to = max(data2$range), by = .25), Mass_.g. = mean(data2$Mass_.g.[data2$Morphotype=="Grey-headed"]), Morphotype = "Grey-headed")
df_PS = data.frame(range = seq(from = min(data2$range), to = max(data2$range), by = .25), Mass_.g. = mean(data2$Mass_.g.[data2$Morphotype=="Pink-sided"]), Morphotype = "Pink-sided")
df_OR= data.frame(range = seq(from = min(data2$range), to = max(data2$range), by = .25), Mass_.g. = mean(data2$Mass_.g.[data2$Morphotype=="Oregon"]), Morphotype = "Oregon")
df_SC= data.frame(range = seq(from = min(data2$range), to = max(data2$range), by = .25), Mass_.g. = mean(data2$Mass_.g.[data2$Morphotype=="Slate-colored"]), Morphotype = "Slate-colored")
df_YE = data.frame(range = seq(from = min(data2$range), to = max(data2$range), by = .25), Mass_.g. = mean(data2$Mass_.g.[data2$Morphotype=="Yellow-eyed"]), Morphotype = "Yellow-eyed")

df_GH$prediction<-predict(lm(Msum_.ml.O2.min.~Mass_.g. + Morphotype * range, data2), newdata=df_GH, type="response", se.fit=F )
df_PS$prediction<-predict(lm(Msum_.ml.O2.min.~Mass_.g. + Morphotype * range, data2), newdata=df_PS, type="response", se.fit=F )
df_OR$prediction<-predict(lm(Msum_.ml.O2.min.~Mass_.g. + Morphotype * range, data2), newdata=df_OR, type="response", se.fit=F )
df_SC$prediction<-predict(lm(Msum_.ml.O2.min.~Mass_.g. + Morphotype * range, data2), newdata=df_SC, type="response", se.fit=F )
df_YE$prediction<-predict(lm(Msum_.ml.O2.min.~Mass_.g. + Morphotype * range, data2), newdata=df_YE, type="response", se.fit=F )

#Plot each morphotype separately
pal <- wes_palette(5, name = "Rushmore1", type = "continuous")
par(mfrow=c(2,3))
plot(Msum_.ml.O2.min.~range, data2[data2$Morphotype=="Grey-headed",], pch=19, col=pal[1],xlim=c(7,23), ylim=c(1,11), xlab="Temperature Annual Range", ylab="Msum")
abline(lm(prediction~range, df_GH), col=pal[1], lwd=2)
title("Grey-headed")

plot(Msum_.ml.O2.min.~range, data2[data2$Morphotype=="Pink-sided",], pch=19, col=pal[4],xlim=c(7,23), ylim=c(1,11), xlab="Temperature Annual Range", ylab="Msum")
abline(lm(prediction~range, df_PS), col=pal[4], lwd=2)
title("Pink-sided")

plot(Msum_.ml.O2.min.~range, data2[data2$Morphotype=="Oregon",], pch=19, col=pal[2],xlim=c(7,23), ylim=c(1,11), xlab="Temperature Annual Range", ylab="Msum")
abline(lm(prediction~range, df_OR), col=pal[2], lwd=2)
title("Oregon")

plot(Msum_.ml.O2.min.~range, data2[data2$Morphotype=="Slate-colored",], pch=19, col=pal[5],xlim=c(7,23), ylim=c(1,11), xlab="Temperature Annual Range", ylab="Msum")
abline(lm(prediction~range, df_SC), col=pal[5], lwd=2)
title("Slate-colored")

plot(Msum_.ml.O2.min.~range, data2[data2$Morphotype=="Yellow-eyed",], pch=19, col=pal[3],xlim=c(7,23), ylim=c(1,11), xlab="Temperature Annual Range", ylab="Msum")
abline(lm(prediction~range, df_YE), col=pal[3], lwd=2)
title("Yellow-eyed")
