#**********************************************************************************************************************************
#**********************************************************************************************************************************

# Project: Pacific White-fronted Goose Population Analysis
# Date: 11 December 2010
# Author: Stephanie Cunningham
# Description: Extracting data from NCEP reanalysis

#**********************************************************************************************************************************
#**********************************************************************************************************************************
rm(list=ls())

library(tidyverse)
library(RNCEP)
library(plotrix)
library(regclass)
library(usdm)
library(psych)

# Read in locations that need a time series
#pts <- read.csv("output/breeding_geocent.csv")
# pts <- read.csv("E:/ResearchProjects/PacificWhite-frontPopulation/output/breeding_geocent.csv")
# pts <- pts[,-1]

# dt <- seq(as.Date("1960-07-01"), as.Date("2017-07-31"), by="day")
# dt <- dt[format(dt, "%m")=="07"]

# dat <- data.frame(date=dt, time="12:00:00", YKD.lat=pts$Latitude[2], YKD.lon=pts$Longitude[2], 
#                   BBL.lat=pts$Latitude[1], BBL.lon=pts$Longitude[1])
# dat <- unite(dat, "datetime", 1:2, sep=" ")


dt <- seq(as.Date("1979-12-01"), as.Date("2017-02-28"), by="day")
dt <- dt[format(dt, "%m")=="12" | format(dt, "%m")=="01" | format(dt, "%m")=="02"]

pts <- data.frame(loc=c("SV"), Latitude=c(39.352642), Longitude=c(-121.978280))

dat <- data.frame(date=dt, time="12:00:00", SV.lat=pts$Latitude[1], SV.lon=pts$Longitude[1])
dat <- unite(dat, "datetime", 1:2, sep=" ")

# temperature
temp <- NCEP.interp(variable="air.sig995", level="surface", lat=dat$SV.lat,
                    lon=dat$SV.lon, dt=dat$datetime, reanalysis2=FALSE, interp="linear")

dat2 <- data.frame(date=dt, SV.temp=temp)
write.csv(dat2, "E:/ResearchProjects/PacificWhite-frontPopulation/output/winter_temp_data.csv")

# Precipitable water
# prwa <- NCEP.interp(variable="pr_wtr.eatm", level="surface", lat=c(dat$YKD.lat, dat$BBL.lat), lon=c(dat$YKD.lon, dat$BBL.lon), dt=dat$datetime, reanalysis2=FALSE, interp="linear")
# dat2 <- data.frame(date=dt, YKD.prwa=prwa[1:1798], BBL.prwa=prwa[1799:3596])
# write.csv(dat2, "E:/ResearchProjects/PacificWhite-frontPopulation/output/prwa_data.csv")

# Relative humidity
# rehu <- NCEP.interp(variable="rhum.sig995", level="surface", lat=c(dat$YKD.lat, dat$BBL.lat), lon=c(dat$YKD.lon, dat$BBL.lon), dt=dat$datetime, reanalysis2=FALSE, interp="linear")
# dat2 <- data.frame(date=dt, YKD.relhum=rehu[1:1798], BBL.relhum=rehu[1799:3596])
# write.csv(dat2, "E:/ResearchProjects/PacificWhite-frontPopulation/output/rehu_data.csv")

# Precipitation rate
prate <- NCEP.interp(variable="prate.sfc", level="gaussian", lat=dat$SV.lat, lon=dat$SV.lon, dt=dat$datetime, reanalysis2=FALSE, interp="linear")
dat2 <- data.frame(date=dt, SV.prate=prate)
write.csv(dat2, "E:/ResearchProjects/PacificWhite-frontPopulation/output/winter_prate_data.csv")

temp <- read.csv("output/winter_temp_data.csv")
# prwa <- read.csv("output/prwa_data.csv")
# rehu <- read.csv("output/rehu_data.csv")
prate <- read.csv("output/winter_prate_data.csv")

# Combine weather by years
dat2 <- data.frame(date=dt, SV.temp=temp$SV.temp, 
                   SV.prate=prate$SV.prate)
dat2$year <- as.numeric(format(dat2$date, "%Y"))
dat2$month <- as.numeric(format(dat2$date, "%m"))
dat2$year[dat2$month==12] <- dat2$year[dat2$month==12] + 1 # Increase the year from the previous December to aggregate winter
dat2 <- dat2[,-5]

sv.avg <- aggregate(cbind(SV.temp, SV.prate)~year, data=dat2, FUN=mean)
sv.ste <- aggregate(cbind(SV.temp, SV.prate)~year, data=dat2, FUN=std.error)

sv.uci.temp <- sv.avg$SV.temp + 1.96*sv.ste$SV.temp
sv.lci.temp <- sv.avg$SV.temp - 1.96*sv.ste$SV.temp
sv.uci.prate <- sv.avg$SV.prate + 1.96*sv.ste$SV.prate
sv.lci.prate <- sv.avg$SV.prate - 1.96*sv.ste$SV.prate

sv.covar <- data.frame(year=as.numeric(sv.avg$year), sv.temp=sv.avg$SV.temp, sv.temp.lcl=sv.lci.temp, sv.temp.ucl=sv.uci.temp,
                        sv.prate=sv.avg$SV.prate, sv.prate.lcl=sv.lci.prate, sv.prate.ucl=sv.uci.prate)

sv.covar <- sv.covar[sv.covar$year>=1980 & sv.covar$year<=2015,]

write.csv(sv.covar, "output/winterSV_covars_ns.csv")


# Combine weather variables by years
dat2 <- data.frame(date=dt, YKD.temp=temp$YKD.temp, BBL.temp=temp$BBL.temp, YKD.prate=prate$YKD.prate, BBL.prate=prate$BBL.prate)
dat2$year <- format(dat2$date, "%Y")

# Convert tempt to C (from K) and precip rate from kg/m^2/s to mm/h
dat2$YKD.temp <- dat2$YKD.temp - 273
dat2$BBL.temp <- dat2$BBL.temp - 273
dat2$YKD.prate <- dat2$YKD.prate * 3600
dat2$BBL.prate <- dat2$BBL.prate * 3600

ykd.avg <- aggregate(cbind(YKD.temp, YKD.prate)~year, data=dat2, FUN=mean)
ykd.ste <- aggregate(cbind(YKD.temp, YKD.prate)~year, data=dat2, FUN=std.error)
ykd.std <- aggregate(cbind(YKD.temp, YKD.prate)~year, data=dat2, FUN=sd)
ykd.uci.temp <- ykd.avg$YKD.temp + 1.96*ykd.ste$YKD.temp
ykd.lci.temp <- ykd.avg$YKD.temp - 1.96*ykd.ste$YKD.temp

ykd.uci.prate <- ykd.avg$YKD.prate + 1.96*ykd.ste$YKD.prate
ykd.lci.prate <- ykd.avg$YKD.prate - 1.96*ykd.ste$YKD.prate

ykd.covar <- data.frame(year=as.numeric(ykd.avg$year), ykd.temp=ykd.avg$YKD.temp, ykd.temp.sd=ykd.std$YKD.temp, ykd.temp.lcl=ykd.lci.temp, ykd.temp.ucl=ykd.uci.temp, 
                        ykd.prate=ykd.avg$YKD.prate, ykd.prate.sd=ykd.std$YKD.prate, ykd.prate.lcl=ykd.lci.prate, ykd.prate.ucl=ykd.uci.prate)

ykd.covar <- ykd.covar[ykd.covar$year>=1980 & ykd.covar$year<=2016,]

plot(ykd.avg$YKD.temp, type="l")
lines(ykd.uci.temp, lty=2)
lines(ykd.lci.temp, lty=2)
lines(ykd.avg$YKD.prate, col="blue")
lines(ykd.uci.prate, lty=2, col="blue")
lines(ykd.lci.prate, lty=2, col="blue")

bbl.avg <- aggregate(cbind(BBL.temp, BBL.prate)~year, data=dat2, FUN=mean)
bbl.ste <- aggregate(cbind(BBL.temp, BBL.prate)~year, data=dat2, FUN=std.error)
bbl.std <- aggregate(cbind(BBL.temp, BBL.prate)~year, data=dat2, FUN=sd)
bbl.uci.temp <- bbl.avg$BBL.temp + 1.96*bbl.ste$BBL.temp
bbl.lci.temp <- bbl.avg$BBL.temp - 1.96*bbl.ste$BBL.temp

bbl.uci.prate <- bbl.avg$BBL.prate + 1.96*bbl.ste$BBL.prate
bbl.lci.prate <- bbl.avg$BBL.prate - 1.96*bbl.ste$BBL.prate

bbl.covar <- data.frame(year=as.numeric(bbl.avg$year), bbl.temp=bbl.avg$BBL.temp, bbl.temp.sd=bbl.std$BBL.temp, bbl.temp.lcl=bbl.lci.temp, bbl.temp.ucl=bbl.uci.temp, 
                        bbl.prate=bbl.avg$BBL.prate, bbl.prate.sd=bbl.std$BBL.prate, bbl.prate.lcl=bbl.lci.prate, bbl.prate.ucl=bbl.uci.prate)

bbl.covar <- bbl.covar[bbl.covar$year>=1980 & bbl.covar$year<=2016,]

# Calculate weighted average depending on banding location
locbyyr <- read.csv("output/bandloc-by-year.csv")
locbyyr <- locbyyr[,-1]


locbyyr <- cbind(locbyyr, ykd.covar, bbl.covar)
locbyyr <- locbyyr[,-6]

locbyyr$prop.BBL[is.na(locbyyr$prop.BBL)] <- 0.5
locbyyr$prop.YKD[is.na(locbyyr$prop.YKD)] <- 0.5

locbyyr <- mutate(locbyyr, temp.wm=(prop.BBL*bbl.temp)+(prop.YKD*ykd.temp))
locbyyr <- mutate(locbyyr, temp.sd=(prop.BBL*bbl.temp.sd)+(prop.YKD*ykd.temp.sd))
locbyyr <- mutate(locbyyr, temp.wlcl=(prop.BBL*bbl.temp.lcl)+(prop.YKD*ykd.temp.lcl))
locbyyr <- mutate(locbyyr, temp.wucl=(prop.BBL*bbl.temp.ucl)+(prop.YKD*ykd.temp.ucl))

locbyyr <- mutate(locbyyr, prate.wm=(prop.BBL*bbl.prate)+(prop.YKD*ykd.prate))
locbyyr <- mutate(locbyyr, prate.sd=(prop.BBL*bbl.prate.sd)+(prop.YKD*ykd.prate.sd))
locbyyr <- mutate(locbyyr, prate.wlcl=(prop.BBL*bbl.prate.lcl)+(prop.YKD*ykd.prate.lcl))
locbyyr <- mutate(locbyyr, prate.wucl=(prop.BBL*bbl.prate.ucl)+(prop.YKD*ykd.prate.ucl))

plot(locbyyr$prate.wm, type="l")
lines(locbyyr$prate.wlcl, lty=2)
lines(locbyyr$prate.wucl, lty=2)
lines(locbyyr$temp.wm, col="blue")
lines(locbyyr$temp.wlcl, col="blue", lty=2)
lines(locbyyr$temp.wucl, col="blue", lty=2)

locbyyr <- locbyyr[,c(1,23:30)]
write.csv(locbyyr, "output/weightedavg_temp_prcp.csv")

# read in landscape metrics data
roost <- read.csv("data/roost_metrics.csv", stringsAsFactors=FALSE)

roost$ENN_MNsc <- scale(roost$ENN_MN)
roost$ENN_AMsc <- scale(roost$ENN_AM)
roost$ENN_SDsc <- scale(roost$ENN_SD)

roost <- mutate(roost, ENN_SEsc=ENN_SDsc/sqrt(NP))

plot(1980:2017, roost$ENN_AMsc, pch=16)
segments(1980:2017, y0=roost$ENN_AMsc-roost$ENN_SEsc, y1=roost$ENN_AMsc+roost$ENN_SEsc)
segments(1980:2017, y0=roost$ENN_AMsc-roost$ENN_SDsc, y1=roost$ENN_AMsc+roost$ENN_SDsc)

# write.csv(roost, "output/scaled_ENN_AM.csv")

klamath <- read.csv("data/klamath_pasture.csv")
klamath$CAsc <- scale(klamath$CA)
klamath$PLANDsc <- scale(klamath$PLAND)

plot(klamath$PLANDsc)

write.csv(klamath, "output/klamath_pasture.csv")

s <- read.csv("C:/Users/stcunnin/Box Sync/Manuscripts/PacificFlyway_IPM-Ibis/revision/non-hunting_survival_covars.csv")

cor(s[1:37,2:4])
s.mod <- lm(Year~., data=s)
VIF(s.mod)
vif(s[1:37,2:4])

ar <- read.csv("C:/Users/stcunnin/Box Sync/Manuscripts/PacificFlyway_IPM-Ibis/revision/age-ratio_covars.csv")
ar <- mutate(ar, rice.intx1=rice_standarized*no.burn)
ar <- mutate(ar, rice.intx2=rice_standarized*burn.decrease)
ar <- mutate(ar, rice.intx3=rice_standarized*no_burn_dummy)
ar$pop.est <- scaledpop[1:36]
ar$pop.est[34] <- mean(c(ar$pop.est[33], ar$pop.est[35]))
cor(ar[1:36,c(2,9:12)])
ar.mod <- lm(Year~., data=ar)
VIF(ar.mod)
vif(ar[1:36,c(2,4,13:17)])

# pairs.panels(ar[,c(2,4,13:17)], method="pearson", hist.col="gray50",density=TRUE,lm=TRUE,ellipses=FALSE)
pairs.panels(ar[,c(2,4,7,17:18)], method="pearson", hist.col="gray50",density=TRUE,lm=TRUE,ellipses=FALSE)
vif(ar[1:36,c(2,4,7,17:18)])


burn2 <- read.csv("C:/Users/stcunnin/Box Sync/Manuscripts/PacificFlyway_IPM-Ibis/revision/2time_burn.csv")
burn2$pop.est <- scaledpop[1:36]
burn2$pop.est[34] <- mean(c(burn2$pop.est[33], burn2$pop.est[35]))

vif(burn2[1:36,c(2,8,3)])
pairs.panels(burn2[,c(4,5,8,3,6:7)], method="pearson", hist.col="gray50",density=TRUE,lm=TRUE,ellipses=FALSE)

burn3 <- read.csv("C:/Users/stcunnin/Box Sync/Manuscripts/PacificFlyway_IPM-Ibis/revision/3time_burn_3moBEST.csv")
burn3$pop.est <- scaledpop[1:37]
burn3$pop.est[34] <- mean(c(burn3$pop.est[33], burn3$pop.est[35]))
burn3[,3:5] <- scale(burn3[,3:5])

vif(burn3[1:36,c(2,10,3)])
pairs.panels(burn3[,c(4,5,10,3)], method="pearson", hist.col="gray50",density=TRUE,lm=TRUE,ellipses=FALSE)


vif(burn3[1:36,c(2:5,7,8)])


## PCA
# YKD weather
# ykd.pca <- ykd.avg[21:57,2:5]
# rownames(ykd.pca) <- ykd.avg[21:57,1]
# pr.out <- prcomp(ykd.pca, scale=FALSE)
# pr.out
# 
# biplot(pr.out, scale=0, xlim=c(-4,5))
# plot(pr.out$x[,1], type="l")
# lines(pr.out$x[,2], col="blue")
# lines(pr.out$x[,3], col="green")
# lines(pr.out$x[,4], col="red")
# 
# pr.var=pr.out$sdev^2
# pr.var
# pve=pr.var/sum(pr.var)
# pve
# 
# ykd.out <- as.data.frame(pr.out$x)
# 
# ykd.out$pdo <- pdo.s$MAY
# cor(ykd.out[,c(1:2,5)])
# 
# 
# ykd.out <- cbind(ykd.out, ykd.avg[21:57,2:5])
# 
# plot(ykd.out$PC1, type="l")
# lines(ykd.out$pdo, col="red")
# lines(ykd.out$YKD.temp, col="blue")
# lines(ykd.out$YKD.prate, col="green")
# 
# # BBL weather
# bbl.pca <- bbl.avg[21:57,2:5]
# rownames(bbl.pca) <- bbl.avg[21:57,1]
# pr.out <- prcomp(bbl.pca, scale=FALSE)
# pr.out
# 
# biplot(pr.out, scale=0, xlim=c(-4,5))
# 
# plot(pr.out$x[,1], type="l")
# lines(pr.out$x[,2], col="blue")
# lines(pr.out$x[,3], col="green")
# lines(pr.out$x[,4], col="red")
# 
# pr.var=pr.out$sdev^2
# pr.var
# pve=pr.var/sum(pr.var)
# pve
# bbl.out <- as.data.frame(pr.out$x)
# 
# all <- cbind(bbl.out, ykd.out)
# all <- all[,c(1,5)]
# names(all) <- c("BBL.PC1", "YKD.PC1")
# 
# # Calculate weighted average depending on banding location
# locbyyr <- read.csv("output/bandloc-by-year.csv")
# locbyyr <- locbyyr[,-1]
# 
# locbyyr <- cbind(locbyyr, all)
# locbyyr$PC1 <- locbyyr$YKD.PC1
# locbyyr$PC1[locbyyr$prop.BBL==1] <- locbyyr$BBL.PC1[locbyyr$prop.BBL==1]
# 
# 
# locbyyr$prop.BBL[is.na(locbyyr$prop.BBL)] <- 0.5
# locbyyr$prop.YKD[is.na(locbyyr$prop.YKD)] <- 0.5
# 
# locbyyr <- mutate(locbyyr, weight.PC1=(prop.BBL*BBL.PC1)+(prop.YKD*YKD.PC1))
# plot(locbyyr$weight.PC1, type="l")
# lines(locbyyr$BBL.PC1, col="red")
# lines(locbyyr$YKD.PC1, col="blue")
# 
# locbyyr <- mutate(locbyyr, weightsd.pc1=sqrt(sum(c(prop.BBL, prop.YKD))*(weight.PC1-)   ))
# 
# 
# write.csv(locbyyr, "output/weighted_pca.csv")
# 
# 
