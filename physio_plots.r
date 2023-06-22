library(investr) #for predFit

mo2 = read.csv("C:\\Users\\nickh\\Documents\\OSU\\project\\data\\physiology\\mo2.csv", header=TRUE, as.is=TRUE)
recovery = read.csv("C:\\Users\\nickh\\Documents\\OSU\\project\\data\\physiology\\recovery.csv", header=TRUE, as.is=TRUE)
ctmax = read.csv("C:\\Users\\nickh\\Documents\\OSU\\project\\data\\physiology\\ctmax.csv", header=TRUE, as.is=TRUE)
mo2 = mo2[!mo2$Note=="exclude" & !mo2$Type=="RMR", !names(mo2)%in%c("Note","Descrip")]
ctmax = ctmax[ctmax$Basin!="Deschutes_Odell",]

data = c(mo2$MO2, recovery$Minutes/60, ctmax$Temperature)
trout = c(mo2$Basin, recovery$Basin, ctmax$Basin)
temp = c(mo2$TempTrial1, recovery$TempTrial1, rep(999, length(ctmax$Basin)))
stat = c(mo2$Type, rep("REC", length(recovery$Basin)), rep("CTM", length(ctmax$Basin)))
agg = aggregate(data, list(temp=temp, trout=trout, stat=stat), function(x) c(m=mean(x), se=sd(x)/sqrt(length(x))))
agg = as.data.frame(cbind(agg[,-ncol(agg)], agg$x))

smr_mod_keno = lm(data~I(temp^2)+temp, subset=trout=="Spencer"&stat=="SMR")
smr_mod_lake = lm(data~I(temp^2)+temp, subset=trout=="Williamson"&stat=="SMR")
#mmr_mod_keno = lm(data~I(temp^2)+temp, subset=trout=="Spencer"&stat=="MMR")
#mmr_mod_lake = lm(data~I(temp^2)+temp, subset=trout=="Williamson"&stat=="MMR")
mmr_mod_keno = nls(data~a*exp(-.5*((temp+b)/c)^2), subset=trout=="Spencer"&stat=="MMR", start=list(a=20, b=-18, c=18), lower=c(0, -30, 0), algorithm="port")
mmr_mod_lake = nls(data~a*exp(-.5*((temp+b)/c)^2), subset=trout=="Williamson"&stat=="MMR", start=list(a=20, b=-18, c=18), lower=c(0, -30, 0), algorithm="port")
#aas_mod_keno = lm(data~I(temp^2)+temp, subset=trout=="Spencer"&stat=="AAS")
#aas_mod_lake = lm(data~I(temp^2)+temp, subset=trout=="Williamson"&stat=="AAS")
aas_mod_keno = nls(data~a*exp(-.5*((temp+b)/c)^2), subset=trout=="Spencer"&stat=="AAS", start=list(a=20, b=-18, c=18), lower=c(0, -30, 0), algorithm="port")
aas_mod_lake = nls(data~a*exp(-.5*((temp+b)/c)^2), subset=trout=="Williamson"&stat=="AAS", start=list(a=20, b=-18, c=18), lower=c(0, -30, 0), algorithm="port")
print(summary(smr_mod_keno))
print(summary(smr_mod_lake))
print(summary(mmr_mod_keno))
print(summary(mmr_mod_lake))
print(summary(aas_mod_keno))
print(summary(aas_mod_lake))
new_keno = seq(14, 25, .1)
new_lake = seq(11, 22, .1)
smr_xy_keno = predict(smr_mod_keno, newdata=data.frame(temp=new_keno), interval="confidence")
smr_xy_lake = predict(smr_mod_lake, newdata=data.frame(temp=new_lake), interval="confidence")
mmr_xy_keno = predFit(mmr_mod_keno, newdata=data.frame(temp=new_keno), interval="confidence")
mmr_xy_lake = predFit(mmr_mod_lake, newdata=data.frame(temp=new_lake), interval="confidence")
aas_xy_keno = predFit(aas_mod_keno, newdata=data.frame(temp=new_keno), interval="confidence")
aas_xy_lake = predFit(aas_mod_lake, newdata=data.frame(temp=new_lake), interval="confidence")

#jpeg("ukl_keno_physiology.jpeg", height=6, width=6, res=400, units="in")
tiff("ukl_keno_physiology.tiff", height=6, width=6, res=400, units="in")
par(mfrow=c(2,2), mar=c(4.1, 4.3, 2.1, 1.1))
mmr_smr = stat %in% c("MMR", "SMR")
mmr_smr_agg = agg$stat %in% c("MMR", "SMR")
pcols = c("dodgerblue", "red", adjustcolor(c("dodgerblue", "red"), alpha.f=.2))

plot(NA, NA, #jitter(temp[mmr_smr]), data[mmr_smr], col=pcols[(trout[mmr_smr]=="Spencer")+3], 
	xlim=c(10,26), ylim=range(data[mmr_smr]), pch=24+(stat[mmr_smr]=="SMR"),
	xlab=expression("Temperature ("*degree*"C)"), ylab=expression("Metabolic Rate (g "*O[2]~kg^-1~min^-1*")"))
polygon(c(new_lake, rev(new_lake)), c(smr_xy_lake[,2], rev(smr_xy_lake[,3])), col=pcols[3], border=NA)
polygon(c(new_lake, rev(new_lake)), c(mmr_xy_lake[,2], rev(mmr_xy_lake[,3])), col=pcols[3], border=NA)
polygon(c(new_keno, rev(new_keno)), c(smr_xy_keno[,2], rev(smr_xy_keno[,3])), col=pcols[4], border=NA)
polygon(c(new_keno, rev(new_keno)), c(mmr_xy_keno[,2], rev(mmr_xy_keno[,3])), col=pcols[4], border=NA)
lines(new_lake, smr_xy_lake[,1], col=pcols[1])
lines(new_lake, mmr_xy_lake[,1], col=pcols[1])
lines(new_keno, smr_xy_keno[,1], col=pcols[2])
lines(new_keno, mmr_xy_keno[,1], col=pcols[2])
points(agg$temp[mmr_smr_agg], agg$m[mmr_smr_agg], col=pcols[(agg$trout[mmr_smr_agg]=="Spencer")+1],
	bg=pcols[(agg$trout[mmr_smr_agg]=="Spencer")+1], pch=24+(agg$stat[mmr_smr_agg]=="SMR"))
segments(x0=agg$temp[mmr_smr_agg], y0=(agg$m+agg$se)[mmr_smr_agg], y1=(agg$m-agg$se)[mmr_smr_agg], col=pcols[(agg$trout[mmr_smr_agg]=="Spencer")+1])
mtext("a)", side=3, line=.5, at=10, cex=1.5)

aas = stat %in% c("AAS")
aas_agg = agg$stat %in% c("AAS")
plot(NA, NA, #jitter(temp[aas]), data[aas], col=pcols[(trout[aas]=="Spencer")+3],
	xlim=c(10,26), ylim=range(data[aas]), xlab=expression("Temperature ("*degree*"C)"), ylab=expression("AAS (g "*O[2]~kg^-1~min^-1*")"))
polygon(c(new_lake, rev(new_lake)), c(aas_xy_lake[,2], rev(aas_xy_lake[,3])), col=pcols[3], border=NA)
polygon(c(new_keno, rev(new_keno)), c(aas_xy_keno[,2], rev(aas_xy_keno[,3])), col=pcols[4], border=NA)
lines(new_lake, aas_xy_lake[,1], col=pcols[1])
lines(new_keno, aas_xy_keno[,1], col=pcols[2])
points(agg$temp[aas_agg], agg$m[aas_agg], pch=21, col=pcols[(agg$trout[aas_agg]=="Spencer")+1], bg=pcols[(agg$trout[aas_agg]=="Spencer")+1])
segments(x0=agg$temp[aas_agg], y0=(agg$m+agg$se)[aas_agg], y1=(agg$m-agg$se)[aas_agg], col=pcols[(agg$trout[aas_agg]=="Spencer")+1])
mtext("b)", side=3, line=.5, at=10, cex=1.5)

fas = stat %in% c("FAS")
fas_agg = agg$stat %in% c("FAS")
plot(jitter(temp[fas]), data[fas], col=pcols[(trout[fas]=="Spencer")+3], xlim=c(10,26),
	xlab=expression("Temperature ("*degree*"C)"), ylab="FAS")
abline(h=2:3, col="grey", lty=1:2, lwd=1.5)
points(agg$temp[fas_agg], agg$m[fas_agg], pch=21, col=pcols[(agg$trout[fas_agg]=="Spencer")+1], bg=pcols[(agg$trout[fas_agg]=="Spencer")+1])
segments(x0=agg$temp[fas_agg], y0=(agg$m+agg$se)[fas_agg], y1=(agg$m-agg$se)[fas_agg], col=pcols[(agg$trout[fas_agg]=="Spencer")+1])
mtext("c)", side=3, line=.5, at=10, cex=1.5)

rec = stat %in% c("REC")
rec_agg = agg$stat %in% c("REC")
plot(jitter(temp[rec]), data[rec], col=pcols[(trout[rec]=="Spencer")+3], xlim=c(10,26),
	xlab=expression("Temperature ("*degree*"C)"), ylab="Recovery Time (h)")
points(agg$temp[rec_agg], agg$m[rec_agg], pch=21, col=pcols[(agg$trout[rec_agg]=="Spencer")+1], bg=pcols[(agg$trout[rec_agg]=="Spencer")+1])
segments(x0=agg$temp[rec_agg], y0=(agg$m+agg$se)[rec_agg], y1=(agg$m-agg$se)[rec_agg], col=pcols[(agg$trout[rec_agg]=="Spencer")+1])
legend("topleft", pch=16, col=c("dodgerblue", "red"), border=NA, legend=c("UKL", "Keno"), inset=c(.1, .1))
mtext("d)", side=3, line=.5, at=10, cex=1.5)

dev.off()

#add stat tests for reviewer 2
#two-sided, spencer == 22 : MMR ns, AAS .11, FAS <.01, REC < .01
#two-sided, spencer >= 22 : MMR .12, AAS <.05, FAS <.01, REC <.01
#one-sided, spencer == 22 : MMR ns, SMR = .07, AAS .055, FAS <.01, REC <.01
#one-sided, spencer >= 22 : MMR .06, SMR ns, AAS <.01, FAS <.01, REC <.01
#fit quadratic lm to SMR, three-parameter gaussian nls to MMR, citing Chen 2015 J Exp Biol
#remove background points? and overlay curves
statsIndex = c("SMR", "MMR", "AAS", "FAS", "REC", "CTM")
for(i in 1:length(statsIndex)){
print(statsIndex[i])
print(wilcox.test(data[stat==statsIndex[i]&trout=="Spencer"&temp>=22], #temp for CTM = 999 so all >=22
	data[stat==statsIndex[i]&trout=="Williamson"&temp>=22], 
	alternative=ifelse(statsIndex[i]%in%c("SMR", "REC"), "less", "greater")))
}