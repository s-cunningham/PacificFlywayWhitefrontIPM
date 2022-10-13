#**********************************************************************************************************************************
#**********************************************************************************************************************************

# Project: Pacific White-fronted Goose Population Analysis
# Date: 5 July 2018
# Author: Stephanie Cunningham
# Description: Extracting dataset from GameBirds data records - mark-recovery model

#**********************************************************************************************************************************
#**********************************************************************************************************************************

library(Hmisc)
library(tidyverse)
library(rgdal)
library(rgeos)

rm(list=ls())

# North America shape file
na <- readOGR("data/spatialdata", "ne_50m_land")
na <- spTransform(na, CRS("+proj=longlat +datum=WGS84"))
par(mar=c(4.5,4.5,1,1))
plot(na, xlim=c(-172,-100),ylim=c(30, 70), axes=TRUE, ylab="Latitude", xlab="Longitude")

# Read in breeding grounds shape file
ak <- readOGR("data/spatialdata", "breeding2") 
ak <- spTransform(ak, CRS("+proj=longlat +datum=WGS84"))
plot(ak, add=TRUE, col="green")

# Read in banding data - just summaries of how many birds banded per year
setwd("data/BandingRecords")
bands <- list.files(path="./", pattern='.csv', all.files=TRUE, full.names=FALSE)
bands <- lapply(bands, FUN=read.csv, header=TRUE, stringsAsFactors=FALSE)
bands <- do.call("rbind", bands)
setwd("E:/ResearchProjects/PacificWhite-frontPopulation")

bands <- subset(bands, B.Year>=1980) # Looks like poor-quality data before 1965

# Subset by birds banded on YK Delta (exclude North Slope and interior AK birds - assumed to be Central Flyway Breeders)
bands <- subset(bands, (GISBLat <= 50) | (GISBLat <= 61 & GISBLong <= -153) | (GISBLat <= 63.2 & GISBLong <= -161))
bands <- subset(bands, Region..State=="Alaska")

# Save coordinates of banding locations and plot
xy <- bands[,c(25,24)]
names(xy) <- c("Longitude", "Latitude")
points(xy, col="blue", pch=16)
abline(h=60)
abline(h=58)

xy$loc <- ifelse(xy$Latitude>60, "YKD", "BBL")
xy$loc <- ifelse(xy$Latitude<58, NA, xy$loc)

geocent <- aggregate(cbind(Longitude, Latitude)~loc, data=xy, FUN=mean)
points(geocent[,2:3], col="red")
# write.csv(geocent, "output/breeding_geocent.csv")

plot(na, xlim=c(-172,-140),ylim=c(55, 70), axes=TRUE, ylab="Latitude", xlab="Longitude")
plot(ak, add=TRUE, col="green")
points(xy, col="blue", pch=16)
points(geocent[,2:3], col="red", pch=16)
abline(v=-158)

# xy2 <- bands[bands$GISBLat<60,c(25,24)]
# names(xy2) <- c("Longitude", "Latitude")
# points(xy2, col="red", pch=15)

# Count contributions by bander
banders <- bands %>% group_by(B.Year, Permits..Permittee) %>% summarise(Sum=sum(Count.of.Birds))
# as.data.frame(banders)

# What happens if I run model with justYKD data
bands <- bands[bands$GISBLat > 60,]



# clear out unneeded columns and rename
bands <- bands[,c(1,5,6,8,14,16,17,24,25,33,34)]
names(bands) <- Cs(Month,Year,Flyway,State,Age,Sex,Status,Latitude,Longitude,BandPrefix,NumberBirds)

# order by band year
N.order <- order(bands$Year, decreasing=FALSE)
bands <- bands[N.order,]

# Subset to only keep birds with know  age
bands <- subset(bands, Age!="Unknown")
# bands <- subset(bands, Sex!="Unknown")
bands$Sex[bands$Sex=="Female; from subsequent encounter"] <- "Female"

# Reduce age groups to HY and AHY only
bands$Age[bands$Age=="Local"] <- "Hatch Year"
bands$Age[bands$Age=="After Second Year" | bands$Age=="Second Year"] <- "After Hatch Year"
bands$Age[bands$Age=="Hatch Year"] <- "HY"
bands$Age[bands$Age=="After Hatch Year"] <- "AHY"

# How many total bands?
sum(bands$NumberBirds)
aggregate(NumberBirds ~ Age + Year, data=bands, FUN=sum)
band.count.pf <- aggregate(NumberBirds ~ Year, data=bands, FUN=sum)
sum(band.count.pf$NumberBirds)

points(bands$Longitude, bands$Latitude, col="purple", pch=9)



bands$origin <- NA
bands$origin <- ifelse(bands$Latitude>60, "YKD", "BBL")

locbyyr <- aggregate(NumberBirds ~ Year + origin, data=bands, FUN=sum)
locbyyr <- spread(locbyyr, key="origin", value=NumberBirds)

yrs <- data.frame(Year=seq(1980, 2016, 1))

locbyyr <- left_join(yrs, locbyyr, by="Year")
locbyyr$BBL[is.na(locbyyr$BBL)] <- 0
locbyyr$YKD[is.na(locbyyr$YKD)] <- 0

locbyyr <- mutate(locbyyr, prop.BBL=BBL/(BBL+YKD))
locbyyr <- mutate(locbyyr, prop.YKD=YKD/(BBL+YKD))

# write.csv(locbyyr, "output/bandloc-by-year.csv")


# Save bands data frame
setwd("E:/ResearchProjects/PacificWhite-frontPopulation")
write.csv(bands, "output/20200203_bands.csv")

# Read in recovery data
setwd("data/BandEncounters")
recov <- list.files(path="./", pattern='.csv', all.files=TRUE, full.names=FALSE)
recov <- lapply(recov, FUN=read.csv, header=TRUE, stringsAsFactors=FALSE)
recov <- do.call("rbind",recov)
setwd("E:/ResearchProjects/PacificWhite-frontPopulation")

recov <- subset(recov, How.Obt..VHow=="Shot.")

# Remove any duplicates
recov <- unique(recov)

# Clear out unneeded columns
recov <- recov[,c(1,5,7,9,15,16,20,22,25,36,43:46,61,67,69)]
names(recov) <- Cs(Band,OriginalNewBand,BandMonth,RecoveryMonth,BandYear,RecoveryYear,BandFlyway,RecoveryFlyway,BandState,RecoveryState,BandLat,BandLong,RecoveryLat,RecoveryLong,AgeAtBand,Sex,RecoveryCondition)

# Filter dataset
recov <- subset(recov, OriginalNewBand=="Original Band")
# recov <- subset(recov, Sex!="Unknown") # loss of 1097
recov <- subset(recov, AgeAtBand!="Unknown") # loss of 335
recov <- subset(recov, (BandYear>=1980 & BandYear<=2016))
recov <- subset(recov, RecoveryCondition=="Dead")
recov <- subset(recov, RecoveryState!="Saskatchewan")

# Subset by birds banded on YK Delta/ Bristol Bay lowlands (exclude North Slope and interior AK birds - assumed to be Central Flyway Breeders)
recov <- subset(recov, (BandLat <= 50) | (BandLat <= 61 & BandLong <= -153) | (BandLat <= 63.2 & BandLong <= -161))
recov <- subset(recov, BandFlyway=="Alaska")

# Subset recoveries to fall/winter/spring (September - March)
recov <- subset(recov, (RecoveryLong <= -103 & RecoveryLat <= 50))
recov <- subset(recov, RecoveryMonth>=09 | RecoveryMonth<=06)

# What happens of we remove birds banded at BBL
recov <- recov[recov$BandLat > 60,]


# order by band year
N.order <- order(recov$BandYear, decreasing=FALSE)
recov <- recov[N.order,]

# Check birds with no recovery state listed
# nostate <- recov[recov$RecoveryState=="",]
# xy.ns <- nostate[,c(14,13)]
# names(xy.ns) <- c("Longitude", "Latitude")
# points(xy.ns, col="red", pch=1)

# Save coordinates of recovery locations and plot
plot(na, xlim=c(-172,-100),ylim=c(20, 70), axes=TRUE, ylab="Latitude", xlab="Longitude")
xy <- recov[,c(14,13)] # Recover locations
names(xy) <- c("Longitude", "Latitude")
points(xy, col="purple", pch=1)

xy2 <- recov[,c(12,11)] # Banding locations of recovered individuals
names(xy2) <- c("Longitude", "Latitude")
points(xy2, col="red", pch=0)
legend(-180, 20, c("Banding Locations", "Recovery Locations", "Locations where recovered birds were banded"), col=c("blue", "red", "green"), pch=c(16,1,0))

rarea <- data.frame(lat=c(39.352642, 28.682925), lon=c(-121.978280, -107.007600))
points(rarea$lon, rarea$lat, col="green", pch=8)

recov <- unique(recov)

setwd("E:/ResearchProjects/PacificWhite-frontPopulation")
write.csv(recov, "output/20200203_recoveries.csv")

# Assign number to year
yrs <- data.frame(year=seq(1980,2016,1), occasion=seq(1,37,1))

recov$BandOcc <- rep(NA, nrow(recov))
recov$RecovOcc <- rep(NA, nrow(recov))
for (i in 1:nrow(recov)) {
  recov$BandOcc[i] <- recov$BandYear[i]-1980+1  
  recov$RecovOcc[i] <- recov$RecoveryYear[i]-1980+1
}

recov <- subset(recov, RecoveryYear<=2016)

allyrs <- unique(c(recov$BandYear, recov$RecoveryYear))

# Condense age groups into HY and AHY
recov$AgeAtBand[recov$AgeAtBand=="After Hatch Year" | recov$AgeAtBand=="Second Year" | recov$AgeAtBand=="After Second Year"] <- "AHY"
recov$AgeAtBand[recov$AgeAtBand=="Hatch Year" | recov$AgeAtBand=="Local"] <- "HY"
recov$Sex[recov$Sex=="Female; from subsequent encounter"] <- "Female"

# How many birds were recovered?
recov %>% group_by(AgeAtBand) %>% count()
recovyears <- recov %>% group_by(AgeAtBand, RecoveryYear) %>% count()

# save a row indicating when each band was fitted and recovered
brdata <- matrix(NA, ncol=4, nrow=nrow(recov))
for (i in 1:nrow(recov)) {
  brdata[i,1] <- i
  brdata[i,2] <- recov$AgeAtBand[i]
  brdata[i,3] <- recov$BandOcc[i]
  brdata[i,4] <- recov$RecovOcc[i]
}

# Convert to data frame and change data type
brdata <- as.data.frame(brdata)
names(brdata) <- c("Individual", "Age", "BandOcc", "RecovOcc")
brdata <- brdata %>% mutate_all(as.character) 
brdata[c(1,3,4)] <- lapply(brdata[c(1,3,4)], as.numeric)

# Subset by age group
br.ahy <- subset(brdata, Age=="AHY")
br.hy <- subset(brdata, Age=="HY")

# Convert to capture history
ch.ahy <- matrix(0, ncol=nrow(yrs), nrow=nrow(br.ahy))
for (i in 1:nrow(br.ahy)) {
  b <- br.ahy$BandOcc[i]
  r <- br.ahy$RecovOcc[i]
  ch.ahy[i,b] <- 1
  ch.ahy[i,r] <- 1
}

ch.hy <- matrix(0, ncol=nrow(yrs), nrow=nrow(br.hy))
for (i in 1:nrow(br.hy)) {
  b <- br.hy$BandOcc[i]
  r <- br.hy$RecovOcc[i]
  ch.hy[i,b] <- 1
  ch.hy[i,r] <- 1
}

# How many birds banded in each year were recovered?
marked <- bands %>% group_by(Year, Age) %>% summarize(Marked=sum(NumberBirds))
m1 <- bands %>% group_by(Age, Sex) %>% summarize(Marked=sum(NumberBirds))
m2 <- bands %>% group_by(Year, Age) %>% count()

marked$BandOcc <- rep(NA, nrow(marked))
for (i in 1:nrow(marked)) {
  marked$BandOcc[i] <- marked$Year[i]-1980+1
}
marked <- as.data.frame(marked)
marked <- marked[,c(1,4,2,3)]

rdata <- brdata %>% group_by(Age, BandOcc) %>% count()
rdata <- as.data.frame(rdata)
names(rdata)[3] <- "recovered"

N.order <- order(rdata$BandOcc, decreasing=FALSE)
rdata <- rdata[N.order,]

mr <- left_join(marked, rdata, by=c("Age", "BandOcc"))
mr$recovered[is.na(mr$recovered)] <- 0
mr <- mutate(mr, MNR=Marked-recovered) # Subtract the number of recovered birds from the number of marked birds

mnr <- matrix(NA, nrow=sum(mr$MNR), ncol=3)
mnr[,1] <- seq(1,sum(mr$MNR),1)

age <- character()
id <- numeric()
for (i in 1:nrow(mr)) {
  new.age <- rep(mr$Age[i], mr$MNR[i])
  age <- c(age, new.age)
  new.id <- rep(mr$BandOcc[i], mr$MNR[i])
  id <- c(id, new.id)
}

mnr[,2] <- age
mnr[,3] <- id
mnr <- as.data.frame(mnr)
names(mnr) <- c("Individual", "Age", "BandOcc")
mnr[1:3] <- lapply(mnr[1:3], as.character)
mnr[c(1,3)] <- lapply(mnr[c(1,3)], as.numeric)

# Subset by age group
ahy.m <- subset(mnr, Age=="AHY")
hy.m <- subset(mnr, Age=="HY")

# Convert to capture history
ch2.ahy <- matrix(0, ncol=nrow(yrs), nrow=nrow(ahy.m))
for (i in 1:nrow(ahy.m)) {
  b <- ahy.m$BandOcc[i]
  ch2.ahy[i,b] <- 1
}

ch2.hy <- matrix(0, ncol=nrow(yrs), nrow=nrow(hy.m))
for (i in 1:nrow(hy.m)) {
  b <- hy.m$BandOcc[i]
  ch2.hy[i,b] <- 1
}

CH.AHY <- rbind(ch.ahy, ch2.ahy)
CH.HY <- rbind(ch.hy, ch2.hy)

# Find out which years had birds banded
get.first <- function(x) min(which(x!=0))
f.ahy <- apply(CH.AHY, 1, get.first) 
f.hy <- apply(CH.HY, 1, get.first)

# Which years are missing?
yrs <- seq(1,37,1) # Banding occasions
ahy.bo <- unique(f.ahy) # Unique banding occasions present
hy.bo <- unique(f.hy) 

m.ahy <- yrs[!(yrs %in% ahy.bo)] # which banding occasions do not appear in the dataset for AHY birds?
m.hy <- yrs[!(yrs %in% hy.bo)] # which banding occasions do not appear in the dataset for HY birds?

# Create matrix for adding in missing birds (to later remove)
ahy.add <- matrix(0, nrow=length(m.ahy), ncol=37)
for (i in 1:length(m.ahy)) {
  b <- m.ahy[i]
  ahy.add[i,b] <- 1
}

hy.add <- matrix(0, nrow=length(m.hy), ncol=37)
for (i in 1:length(m.hy)) {
  b <- m.hy[i]
  hy.add[i,b] <- 1
}

CH.AHY <- rbind(CH.AHY, ahy.add)
CH.HY <- rbind(CH.HY, hy.add)

f.ahy <- apply(CH.AHY, 1, get.first) 
CH.AHY <- cbind(f.ahy,CH.AHY)
N.order <- order(CH.AHY[,1], decreasing=FALSE)
CH.AHY <- CH.AHY[N.order,]
CH.AHY <- CH.AHY[,2:ncol(CH.AHY)]

f.hy <- apply(CH.HY, 1, get.first)
CH.HY <- cbind(f.hy,CH.HY)
N.order <- order(CH.HY[,1], decreasing=FALSE)
CH.HY <- CH.HY[N.order,]
CH.HY <- CH.HY[,2:ncol(CH.HY)]

# y1 <- matrix(c(1, rep(0, 36)),  nrow=1, ncol=37)

# CH.AHY <- rbind(y1, CH.AHY)
# CH.HY <- rbind(y1, CH.AHY)

# write.csv(CH.HY, "data/ch_hy.csv")
# write.csv(CH.AHY, "data/ch_ahy.csv")

# Define function to create an m-array for mark-recovery (MR) data ... from Kery & Schaub 2012, Chapter 8
marray.dead <- function(MR){
  nind <- dim(MR)[1]
  n.occasions <- dim(MR)[2]
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
  # Create vector with occasion of marking 
  get.first <- function(x) min(which(x!=0))
  f <- apply(MR, 1, get.first)
  # Calculate the number of released individuals at each time period
  first <- as.numeric(table(f))
  for (t in 1:n.occasions){
    m.array[t,1] <- first[t]
  }
  # Fill m-array with recovered individuals
  rec.ind <- which(apply(MR, 1, sum)==2)
  rec <- numeric()
  for (i in 1:length(rec.ind)){
    d <- which(MR[rec.ind[i],(f[rec.ind[i]]+1):n.occasions]==1)
    rec[i] <- d + f[rec.ind[i]]
    m.array[f[rec.ind[i]],rec[i]] <- m.array[f[rec.ind[i]],rec[i]] + 1
  }
  # Calculate the number of individuals that are never recovered
  for (t in 1:n.occasions){
    m.array[t,n.occasions+1] <- m.array[t,1]-sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}

marr.ahy <- marray.dead(CH.AHY)   # m-array for after hatch year birds
marr.hy <- marray.dead(CH.HY)   # m-array for hatch year birds

all0 <- matrix(0, nrow=1, ncol=37)

marr.ahy[m.ahy,] <- all0
marr.hy[m.hy,] <- all0

sum(rowSums(marr.ahy))
sum(rowSums(marr.hy))

# write.csv(marr.ahy, "data/marr_ahy_noBBL.csv")
# write.csv(marr.hy, "data/marr_hy_noBBL.csv")


