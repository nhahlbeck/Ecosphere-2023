#read in tracking data and select only valid active detections
trackingDB = read.csv("~\\OSU\\project\\data\\telemetry\\tracking_20181101_with_keno_repeats_labeled.csv", header=TRUE, as.is=TRUE)
trackingDB$date = as.Date(trackingDB$date, "%m/%d/%Y")
pdate = as.POSIXlt(trackingDB$date)

keno_act = trackingDB$detection == "keno"
lake_act = trackingDB$detection %in% c("lake", "trib")
#remove station readings b/c conglomeration of temps, plus broken sensor
broken_sensor = paste(trackingDB$radio_freq, trackingDB$radio_id)=="250 17" & trackingDB$date>as.Date("2017-03-20")
lake_act = lake_act & trackingDB$leader != "STATION" & !broken_sensor

#extract environmental temps from dates of interest--summer 2017
#usgs data
envTable = read.csv("~\\OSU\\project\\data\\temperature\\usgs_envvars_linkdam.csv", as.is=TRUE)
envTable$date = as.Date(envTable$date, format="%m/%d/%Y")
datelims = as.Date(c("2017-05-01", "2017-10-01"))
edate = as.POSIXlt(envTable$date)
envInds = envTable$date>=datelims[1] & envTable$date<=datelims[2]

keno_env = read.csv("C:\\Users\\nickh\\Documents\\OSU\\project\\data\\temperature\\usgs_envvars_keno.csv", header=TRUE, as.is=TRUE)
keno_env$date = as.Date(keno_env$date, "%m/%d/%Y")
kdate = as.POSIXlt(keno_env$date)
keno_envInds = keno_env$date>=datelims[1] & keno_env$date<=datelims[2]

#plot fish temps vs environmental temps for those dates
#jpeg("keno_comparison.jpeg", width=8, height=4, units="in", res=400)
tiff("keno_comparison.tiff", width=8, height=4, units="in", res=400)
par(mfrow=1:2)
plot(jitter(pdate$yday[lake_act & pdate$year == 117 & pdate$mon%in%4:8], 2), trackingDB$temp_C[lake_act & pdate$year == 117 & pdate$mon%in%4:8],
	col="dodgerblue", cex=.75, pch=16, xlab="Date", ylab=expression("Temperature ("*degree*"C)"), ylim=c(6,30), xlim=as.POSIXlt(datelims)$yday, xaxt="n")
points(jitter(pdate$yday[keno_act], 2), trackingDB$temp_C[keno_act], col="red", cex=.75, pch=16)
lines(edate$yday[envInds], envTable$tempMax[envInds], col="dodgerblue")
lines(kdate$yday[keno_envInds], keno_env$tempMax[keno_envInds], col="red")
axis(side=1, at=edate$yday[envInds & edate$mday==1], labels=month.abb[unique(edate$mon[envInds])+1])
mtext("a)", side=3, line=.5, at=as.POSIXlt("2017-05-01")$yday, cex=1.5)
#legend("topleft", lty=c(1,2,NA,NA), pch=c(NA,NA,16,16), col=c("red", "red", "dodgerblue", "red"), legend=c("lake max", "keno max", "UKL fish", "Keno fish"), horiz=TRUE, cex=0.75)
#plot histogram of jul-aug tag temps for each population
lakehist = hist(trackingDB$temp_C[lake_act & pdate$mon %in% 6:7], breaks=0:30, plot=FALSE)
lakehist$counts = 100*lakehist$counts/sum(lakehist$counts)
kenohist = hist(trackingDB$temp_C[keno_act & pdate$mon %in% 6:7], breaks=0:30, plot=FALSE)
kenohist$counts = 100*kenohist$counts/sum(kenohist$counts)
plot(kenohist, col="red", border=NA, xlim=c(5,30), ylab="% of Jul-Aug Detections", xlab=expression("Tag Temperature ("~degree~"C)"), main="")
plot(lakehist, col="dodgerblue", border=NA, add=TRUE)
mtext("b)", side=3, line=.5, at=5, cex=1.5)
legend("topleft", inset=c(.1,.1), legend=c("UKL", "Keno"), fill=c("dodgerblue", "red"), border=NA, cex=.8)
dev.off()


#calculate fish temp relative to ambient
anomaly = rep(NA, length(trackingDB$date))
for(i in 1:length(anomaly))
{
	if(lake_act[i])
	{
		anomaly[i] = trackingDB$temp_C[i] - envTable$tempMean[envTable$date==trackingDB$date[i]]
	}
	if(keno_act[i])
	{
		anomaly[i] = trackingDB$temp_C[i] - keno_env$tempMean[keno_env$date==trackingDB$date[i]]
	}
}
cat("90%ile lake", quantile(anomaly[lake_act & pdate$mon%in%6:7], .9, na.rm=TRUE), "\n")
cat("10%ile keno", quantile(anomaly[keno_act & pdate$mon%in%6:7], .1, na.rm=TRUE), "\n")

#check difference between keno and lake ambient mean and max during summer 2017
cat("90%ile lake/keno max diff: ", quantile(envTable$tempMax[envInds&edate$mon%in%6:7] - keno_env$tempMax[keno_envInds&kdate$mon%in%6:7], .9, na.rm=TRUE), "\n")
cat("90%ile lake/keno mean diff: ", quantile(envTable$tempMean[envInds&edate$mon%in%6:7] - keno_env$tempMean[keno_envInds&kdate$mon%in%6:7], .9, na.rm=TRUE), "\n")
cat("90%ile lake/keno min diff: ", quantile(envTable$tempMin[envInds&edate$mon%in%6:7] - keno_env$tempMin[keno_envInds&kdate$mon%in%6:7], .9, na.rm=TRUE), "\n")

#compare ambient oxygen at both sites
#library(rMR)
#keno_opp = DO.unit.convert(keno_env$DOmin, "mg/L", "pct", "kpa", 101.325, keno_env$tempMin)
#lake_opp = DO.unit.convert(envTable$DOmin, "mg/L", "pct", "kpa", 101.325, envTable$tempMin)
#allyears = sort(unique(c(edate$year, kdate$year)))
#hypoxdays = matrix(ncol=4, nrow=length(allyears), dimnames=list(allyears+1900,
#	c("L acute", "L chron", "K acute", "K chron")))
#for(i in 1:length(allyears))
#{
#	lsummer = edate$year==allyears[i] & edate$mon%in%5:9
#	ksummer = kdate$year==allyears[i] & kdate$mon%in%5:9
#	hypoxdays[i,1] = sum(lake_opp[lsummer]<5, na.rm=TRUE)
#	hypoxdays[i,2] = sum(lake_opp[lsummer]<10, na.rm=TRUE)
#	hypoxdays[i,3] = sum(keno_opp[ksummer]<5, na.rm=TRUE)
#	hypoxdays[i,4] = sum(keno_opp[ksummer]<10, na.rm=TRUE)
#}
#print(hypoxdays)


#histogram of summer O2 mins at both sites
jpeg("keno_ukl_oxygen_hist.jpeg")
lakehist = hist(envTable$DOmin[edate$mon %in% 5:8], breaks=0:15)
lakehist$counts = lakehist$counts/sum(lakehist$counts)
kenohist = hist(keno_env$DOmin[kdate$mon %in% 5:8], breaks=0:15)
kenohist$counts = kenohist$counts/sum(kenohist$counts)
plot(kenohist, col=adjustcolor("red", alpha.f=.5), border=NA, xlim=c(0,15), ylab="% of Days from Jun-Sep", xlab=expression("Min O"[2]~"(mg L"^-1*")"), main="")
plot(lakehist, col=adjustcolor("dodgerblue", alpha.f=.5), border=NA, add=TRUE)
#legend("topleft", legend=c("UKL", "Keno"), fill=adjustcolor(c("dodgerblue", "red"), alpha.f=.5))
dev.off()

#put 2017 in context of other years for reviewer 1
#plot 25, 50, and 75 daily quantiles with 2017 overlaid
tqs = aggregate(envTable$tempMax, list(edate$yday), function(x) rep(mean(x, na.rm=TRUE), 3)+c(-1,0,1)*rep(sd(x, na.rm=TRUE), 3))$x
	#function(x) quantile(x, c(0, .5, 1), na.rm=TRUE))$x
jpeg("2017_temp_context.jpeg", width=4, height=4, units="in", res=400)
plot(NA,NA, xlim=c(0, 365), ylim=c(1, 26), xlab="", ylab=expression("Temperature ("*degree*"C)"), xaxt="n")
polygon(c(0:365, 365:0), c(tqs[,1], rev(tqs[,3])), border=NA, col=grey(.7))
lines(0:365, tqs[,2])
lines(0:364, envTable$tempMax[edate$year==117], col="red")
axis(side=1, at=edate$yday[edate$year==117 & edate$mday==1], labels=month.abb, las=3)
dev.off()
#plot deviation of daily max temp from mean relative to sd
#jpeg("2017_temp_anomaly.jpeg", width=4, height=4, units="in", res=400)
#tms = aggregate(envTable$tempMax, list(edate$yday), function(x) c(mean(x, na.rm=TRUE), sd(x, na.rm=TRUE)))$x
#tanom = envTable$tempMax[edate$year==117] - tms[-nrow(tms),1]
#plot(NA, NA, xlim=c(0, 364), ylim=range(tanom), type="l", xaxt="n", xlab="", ylab=expression("Temperature anomaly ("*degree*"C)"))
#polygon(c(0:364, 364:0), c(-tms[-nrow(tms),2], rev(tms[-nrow(tms),2])), col=grey(.7), border=NA)
#abline(h=0, lty=2)
#lines(0:364, tanom)
#axis(side=1, at=edate$yday[edate$year==117 & edate$mday==1], labels=month.abb, las=3)
#dev.off()
