#read in data and use only fish with at least one condition metric
condition_data = read.csv("~\\OSU\\project\\data\\condition\\condition2018.csv", header=TRUE, as.is=TRUE)
condition_data$date = as.Date(condition_data$date, "%m/%d/%Y")
fish = condition_data[(!is.na(condition_data$diet_id) | !is.na(condition_data$bioimp_a)
	| !is.na(condition_data$wt_g)),]
#take out spencer ck--not comparable to adults
fish = fish[!fish$location=="Spencer Ck",]
#take out keno diets but leave in table for condition
fish$diet_id[fish$location=="Keno"] = NA


#create vectors for relevant data sources
date = fish$date
datelist = as.POSIXlt(date)
fish18 = datelist$year==118
fish1718 = datelist$year%in%117:118
loc = fish$location
pop = loc
pop[pop != "Keno"] = "Lake"
pop=factor(pop, levels=c("Lake", "Keno"))
lakefish = pop=="Lake"
loc_specific = tolower(loc)
loc_specific[loc_specific%in%c("ukl", "agency l")] = "lake"
loc_specific[loc_specific=="williamson"] = "lower willy" #no diet samples from above sprague
hab = ifelse(loc_specific == "lake", "lake", "trib")
FL = fish$FL_mm
wt = fish$wt_g
K = 100*1000*wt/FL^3
#for lakefish, late jun (27th in PB) abnormally high
#remove outliers likely errors--two from TNC Jul 27 2018, one from Keno Sept 10 2018, one from Agency May 5 2017
K[K>1.99|K<.5] = NA
bioimp_a = fish$bioimp_a
bioimp_b = fish$bioimp_b
phase_angle = atan(bioimp_b/bioimp_a)*(180/pi)
month = factor(month.abb[datelist$mon+1], levels=month.abb[sort(unique(datelist$mon))+1])
halfmonth = factor(ifelse(datelist$mday<14, "early", "late"))
season = rep(NA, length(loc))
season[!loc%in%c("UKL", "Agency L")&datelist$mon%in%(5:8)&fish1718&lakefish] = "Su"
season[loc%in%c("UKL", "Agency L")&datelist$mon%in%(3:5)&fish1718&lakefish] = "Sp"
season[loc%in%c("UKL", "Agency L")&datelist$mon%in%(7:9)&fish1718&lakefish] = "F"
season = factor(season, levels = c("Sp", "Su", "F"))
spsu = season
spsu[!season%in%c("Sp", "Su")] = NA
spsu = factor(spsu, levels=c("Sp", "Su"))


#group fulton K differently than phase angle and diet--longer timescale of response
kgroup = rep(NA, length(pop))
kgroup[((datelist$mon==4 & datelist$mday>=14) | (datelist$mon==5 & datelist$mday<14)) & lakefish & fish1718] = "Sp"
kgroup[(datelist$mon==8 | (datelist$mon==7 & datelist$mday>=14)) & lakefish & fish1718] = "Su"
kgroup[(datelist$mon==9 & lakefish) & fish1718] = "F"
#take out 2 in lake in september--may have started recovering condition
kgroup[fish18 & datelist$mon==8 & loc=="UKL"] = NA
#pick factor designation to match above
#kgroup = factor(kgroup, levels=c("Begin", "End"))
kgroup = factor(kgroup, levels=c("Sp", "Su", "F"))

#before/after groups for keno and lake comparison (may 14-jun 14 and aug 14-sep 10)
#starting in may to account for variable timing of return from spawning--confounds condition
beginend = rep(NA, length(pop))
beginend[fish18 & (datelist$mon==4 | (datelist$mon==5 & datelist$mday<14))] = "Begin"
beginend[fish18 & (datelist$mon==8 | (datelist$mon==7 & datelist$mday>=14))] = "End"
beginend = factor(beginend)
#take out 2 in lake in september--may have started recovering condition
beginend[fish18 & datelist$mon==8 & loc=="UKL"] = NA

spsuK = kgroup
spsuK[!kgroup%in%c("Sp", "Su")] = NA
spsuK = factor(spsuK, levels=c("Sp", "Su"))

#correct phase angle for FL in 2018--fits data in all years well anyway
pamod = lm(phase_angle[fish18]~FL[fish18])
phase_dev = phase_angle - coef(pamod)[1] - coef(pamod)[2]*FL

#create boxplot for phase angle in lake
jpeg("phase_angle.jpeg", height=2.5, width=4, units="in", res=400)
par(bty="n", mar=c(2.1, 4.5, 1.1, 3.1))#c(2.1,4.5,2.1,3.1))
boxcoords = c(1,3,5)
boxcols = "grey"
bp=boxplot(phase_angle~season, boxwex=1, at=boxcoords, xlab="n",
	col=boxcols, ylab=expression("Phase Angle ("*degree*")"), xaxt="n", cex.lab=1.4, cex.axis=1.4)
#mtext(c("N =",bp$n), side=3, line=.5, at=c(0,boxcoords), cex=1)
mtext(levels(season), line=1, at=boxcoords, cex=1.4, side=1)
dev.off()

#boxplot phase angle including keno before/after summer - standardize to before median
meds = c(median(phase_angle[!lakefish&beginend=="Begin"], na.rm=TRUE),
	median(phase_angle[lakefish&beginend=="Begin"], na.rm=TRUE))
phase_std = 100*(phase_angle - meds[lakefish+1])/meds[lakefish+1]
#K before/after summer - standardize to before median
meds = c(median(K[!lakefish&beginend=="Begin"], na.rm=TRUE),
	median(K[lakefish&beginend=="Begin"], na.rm=TRUE))
Kstd = 100*(K - meds[lakefish+1])/meds[lakefish+1]

#jpeg("phase_angle_and_K_keno.jpeg", width=8, height=4, res=400, units="in")
tiff("phase_angle_and_K_keno.tiff", width=8, height=4, res=400, units="in")
par(mar=c(5,4,2,1), mfrow=1:2)
boxcoords = c(1:2, 4:5)
boxcols = rep(c("dodgerblue", "red"), each=2)
bp=boxplot(Kstd~beginend*pop, boxwex=.8, at=boxcoords, col=boxcols,
	ylab="Condition Factor (% deviation)", xlab="", names=rep(levels(beginend),2))
mtext("a)", side=3, line=.5, at=0, cex=1.5)
mtext(bp$n, side=3, line=.5, at=boxcoords, cex=.8)
bp=boxplot(phase_std~beginend*pop, boxwex=.8, at=boxcoords, col=boxcols,
	ylab="Phase Angle (% deviation)", xlab="", names=rep(levels(beginend), 2))
mtext("b)", side=3, line=.5, at=0, cex=1.5)
mtext(bp$n, side=3, line=.5, at=boxcoords, cex=.8)
legend("bottomright", fill=c("dodgerblue", "red"), legend=c("UKL", "Keno"), inset=c(.1, .1))
dev.off()


jpeg("fultons_k.jpeg", height=2.5, width=4, units="in", res=400)
par(bty="n", mar=c(2.1, 5, 1.1, 3.1))#c(2.1,4.5,2.1,3.1))
boxcoords = c(1,3,5)
boxcols = "grey"
bp=boxplot(K~kgroup, boxwex=1, at=boxcoords, xlab="n",
	col=boxcols, ylab=expression("Fulton's K (100 g cm"^{-3}*")"), xaxt="n", cex.lab=1.4, cex.axis=1.4)
#mtext(c("N =",bp$n), side=3, line=.5, at=c(0,boxcoords), cex=1)
mtext(levels(kgroup), line=1, at=boxcoords, cex=1.4, side=1)
dev.off()

#separate months
#jpeg("phase_angle_keno_month.jpeg")
par(mar=c(5,4,2,1))
boxcoords = c(1:12, 15:26)#c(1:6, 9:14)
boxcols = rep(c("grey", "white"), each=12)
bp=boxplot(phase_angle[fish18]~halfmonth[fish18]*droplevels(month[fish18])*pop[fish18], boxwex=.8, at=boxcoords, las=3,
	col=boxcols, ylab="Phase Angle (deg)", xlab="")#, names=rep(levels(droplevels(month[fish18])),2))
mtext(c("N =",bp$n), side=3, line=.5, at=c(0,boxcoords), cex=.8)
legend("topright", fill=c("grey", "white"), legend=c("Lake", "Keno"))
dev.off()

#jpeg("fultons_k_keno_month.jpeg")
par(mar=c(5,4,2,1))
bp=boxplot(K[fish18]~halfmonth[fish18]*droplevels(month[fish18])*pop[fish18], boxwex=.8, at=boxcoords, col=boxcols, las=3,
	ylab=expression("Fulton's K (100 * g/cm"^{3}~")"), xlab="")#, names=rep(levels(droplevels(month[fish18])),2))
mtext(c("N =",bp$n), side=3, line=.5, at=c(0,boxcoords), cex=.8)
legend("topright", fill=c("grey", "white"), legend=c("Lake", "Keno"))
dev.off()

#test significance of before/after summer
print("Phase Angle: UKL")
print(wilcox.test(phase_angle[lakefish&beginend=="Begin"], phase_angle[lakefish&beginend=="End"]))
print("Phase Angle: Keno")
print(wilcox.test(phase_angle[!lakefish&beginend=="Begin"], phase_angle[!lakefish&beginend=="End"]))
print("Fulton's K: UKL")
print(wilcox.test(K[lakefish&beginend=="Begin"], K[lakefish&beginend=="End"]))
print("Fulton's K: Keno")
print(wilcox.test(K[!lakefish&beginend=="Begin"], K[!lakefish&beginend=="End"]))

#stat tests for lake difference
#cat("----------\nSTAT TESTS\n----------\nPA: \n")
#print(pairwise.wilcox.test(phase_angle[lakefish], season[lakefish], "BH"))
#cat("\nK: \n")
#print(pairwise.wilcox.test(K[lakefish&fish1718], kgroup[lakefish&fish1718], "BH"))

#print medians
#cat("--------------------\nKENO DESCRIPTIVE COMPARISON\n--------------------\nPA stats:\n")
#print(by(phase_std[beginend=="End"], list(pop[beginend=="End"]), median, na.rm=TRUE))
#cat("\nK stats:\n")
#print(by(Kstd[beginend=="End"], list(pop[beginend=="End"]), median, na.rm=TRUE))

#confidence intervals for lake values
cint = function(x, a=.95)
{
	cnames = paste0(c(100*a, "", 100*a), c("% lower", "mean", "% upper"))
	a = ifelse(a > .5, 1-a, a)
	m = mean(x, na.rm=TRUE)
	clims = m + qnorm(c(a/2, 1-a/2)) * sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x)))
	cint = c(clims[1], m, clims[2])
	names(cint) = cnames
	return(cint)
}
cat("\nPA CIs:\n")
print(by(phase_angle, season, function(x) round(cint(x, .95), 1)))
cat("\nK CIs:\n")
print(by(K, kgroup, function(x) round(cint(x, .95), 2)))