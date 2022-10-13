
library(tidyverse)
library(regclass)
library(usdm)
library(psych)


ak.weather <- read.csv("output/weightedavg_temp_prcp.csv")
ak.weather <- ak.weather[,-1]
ak.weather <- ak.weather[ak.weather$Year<=2015,]

rice <- read.csv("data/NASSrice.csv")
rice <- rice[1:36,c(2,11)] 

best.raw <- read.csv("data/BEST_data.csv")
best <- data.frame(year=best.raw$YEAR, Dec=c(NA,best.raw$DEC[1:70]), Jan=best.raw$JAN, Feb=best.raw$FEB, DJF.mean=best.raw$Wint.DJF.mean)
best <- best[best$year>=1980 & best$year<=2015,]
best.sd <- best %>% group_by(year) %>% summarise(DJF.sd=sd(c(Dec,Jan,Feb)))
best <- cbind(best, best.sd$DJF.sd)
names(best)[6] <- "DJF.sd"
best <- best[1:36,]

# Plot covariates
par(mfrow=c(2,2), mar=c(2,4.5,1,0.5), cex.axis=1.1, cex.lab=1.5, las=1)
plot(rice$Year, rice$planted.ha/1000, ylim=c(50, 300), type="l", frame.plot=TRUE, lwd=2, cex.lab=1.5, ylab="Rice Acreage (1000 ha)", xlab="")
points(rice$Year[1:20], rice$planted.ha[1:20]/1000, pch=1, cex=2)
points(rice$Year[21:36], rice$planted.ha[21:36]/1000, pch=16, cex=2)
axis(1,labels=F)
text(x=1980, y=300, "A", cex=1.5)
legend("bottomright", c("without burning restrictions", "with burning restrictions"), pch=c(1,16), lty=1, lwd=2, pt.cex=2, bty="n", cex=1.2)

plot(best$year, best$DJF.mean, type="n", ylim=c(-3,3),  cex.lab=1.5, ylab="BEST Index", xlab="")
axis(1,labels=F)
# polygon(x=c(1980:2015,2015:1980), y=c(best$DJF.mean-best$DJF.sd, rev(best$DJF.mean+best$DJF.sd)), col="#cccccc", border="#cccccc")
lines(best$year, best$DJF.mean,lwd=2)
text(x=1980, y=3, "B", cex=1.5)

par(mar=c(2.5,4.5,1,0.5))
plot(ak.weather$Year, ak.weather$temp.wm, type="n", ylim=c(6,14), cex.lab=1.5, xlab="", ylab="July Temperature (C)")
# polygon(x=c(1980:2015,2015:1980), y=c(ak.weather$temp.wm-ak.weather$temp.sd, rev(ak.weather$temp.wm+ak.weather$temp.sd)), col="#cccccc", border="#cccccc")
lines(ak.weather$Year, ak.weather$temp.wm,lwd=2)
text(x=1980, y=14, "C", cex=1.5)

plot(ak.weather$Year, ak.weather$prate.wm, type="l", ylim=c(0,0.15), lwd=2, cex.lab=1.5, xlab="", ylab="July Precipitation Rate (mm/h)")
# polygon(x=c(1980:2015,2015:1980), y=c(ak.weather$prate.wm-ak.weather$prate.sd, rev(ak.weather$prate.wm+ak.weather$prate.sd)), col="#cccccc", border="#cccccc")
lines(ak.weather$Year, ak.weather$prate.wm,lwd=2)
text(x=1980, y=0.15, "D", cex=1.5)
box(which="plot", lty="solid")
title(xlab='Year', cex.lab=2.2, line=1.2, outer=T)










