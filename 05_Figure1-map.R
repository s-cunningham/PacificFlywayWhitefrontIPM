#########################################################################################################################
##  Code for making Fig. 1 in Cunningham et al. 2020. Ibis. 

#########################################################################################################################

library(adehabitatLT)
library(sp)
library(rgdal)
library(tidyverse)
library(spatialEco)
library(ggspatial)
library(ggmap)
library(gridExtra)


rm(list=ls())

# Import polygon of north america
world <- readOGR("E:/ResearchProjects/GPStracks/data/shapefiles", "ne_50m_land") 
world <- spTransform(world, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
world <- spTransform(world, CRS("+proj=aeqd +lat_0=52 +lon_0=-100"))  # North America
# world <- spTransform(world, CRS("+proj=laea +lat_0=45 +lon_0=-100"))  # North America
mp <- borders(world, colour="gray50", fill="gray80")
# ggplot() + mp + coord_cartesian(xlim=c(-2500000,3500000), ylim=c(-3000000,3000000)) + theme_bw() 

outline <- borders(world, colour="gray50", fill=NA)

states <- readOGR("data/spatialdata", "ne_50m_admin_1_states_provinces")
states <- spTransform(states, CRS("+proj=aeqd +lat_0=52 +lon_0=-100"))  
statesb <- borders(states, colour="gray70", fill="gray80")

caneor10m <- readOGR("data/spatialdata", "CANEOR_10m")
caneor10m <- spTransform(caneor10m, CRS("+proj=aeqd +lat_0=52 +lon_0=-100"))  
caneor10mb <- borders(caneor10m, colour="gray50", fill="gray80")

ca10m <- readOGR("data/spatialdata", "California_10m")
ca10m <- spTransform(ca10m, CRS("+proj=aeqd +lat_0=52 +lon_0=-100"))  
ca10mb <- borders(ca10m, colour="gray50", fill="gray90")

sv <- readOGR("data/spatialdata/SacramentoValley", "SacramentoValley1")
sv <- remove.holes(sv)
sv <- spTransform(sv, CRS("+proj=aeqd +lat_0=52 +lon_0=-100")) 
svb <- borders(sv, colour="black", fill="gray80")

rice <- readOGR("data/spatialdata/Projected_AlbersEqualArea", "SacValleyRice")
rice <- spTransform(rice, CRS("+proj=aeqd +lat_0=52 +lon_0=-100"))  
riceb <- borders(rice, colour="gray30", fill="gray30")

bands <- read.csv("output/20200203_bands.csv")
bands <- bands[,-1]

recov <- read.csv("output/20200203_recoveries.csv")
recov <- recov[,-1]

latlon <- data.frame(latitude=c(bands$Latitude, recov$RecoveryLat), longitude=c(bands$Longitude, recov$RecoveryLong))
latlon <- cbind(rep(1,nrow(latlon)), latlon)
names(latlon)[1] <- "group"
aggregate(cbind(latitude,longitude)~group, data=latlon, FUN=mean)


xy <- bands[,c(9,8)]
bands.xy <- SpatialPointsDataFrame(coords=xy, data=bands, proj4string=CRS("+proj=longlat +datum=WGS84"))

# writeOGR(obj=bands.xy, dsn="data/spatialdata", layer="bandsxy", driver="ESRI Shapefile")


bands.xy <- spTransform(bands.xy, CRS("+proj=aeqd +lat_0=52 +lon_0=-100"))
bands.coords <- as.data.frame(bands.xy@coords)
bands.coords <- cbind(rep("Ringing Locations", nrow(bands.coords)), bands.coords)
names(bands.coords) <- c("group", "x", "y")

xy <- recov[,c(14,13)]
recov.xy <- SpatialPointsDataFrame(coords=xy, data=recov, proj4string=CRS("+proj=longlat +datum=WGS84"))
# writeOGR(obj=recov.xy, dsn="data/spatialdata", layer="recovxy", driver="ESRI Shapefile")

recov.xy <- spTransform(recov.xy, CRS("+proj=aeqd +lat_0=52 +lon_0=-100"))
recov.coords <- as.data.frame(recov.xy@coords)
recov.coords <- cbind(rep("Recovery Locations", nrow(recov.coords)), recov.coords)
names(recov.coords) <- c("group", "x", "y")


# Panel A
coords <- rbind(bands.coords, recov.coords)

panelA <- ggplot() + mp + statesb + outline + coord_cartesian(xlim=c(-4000000,0), ylim=c(-2000000,3000000)) + theme_classic() +
  geom_point(data=coords, aes(x=x, y=y, shape=group), size=4) +  
  scale_shape_manual(values=c(0,1))+
  theme(axis.title.x=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
        axis.title.y=element_blank(), panel.border=element_rect(color="black", fill=NA, size=0.5),
        legend.justification=c(0,0), legend.position=c(0,0.5), legend.text=element_text(size=12), 
        legend.title=element_blank(),
        legend.background=element_rect(fill=NA, colour=NA)) +
  annotation_scale(aes(location="bl"), plot_unit="m", text_cex=1) + 
  annotate(geom="text", x=-400000, y=-200000, label="North\nAmerica", fontface="bold", size=6) +
  annotate(geom="text", x=-4000000, y=3100000, label="(a)", fontface="bold", size=6) +
  annotate(geom="text", x=-3400000, y=-400000, label="Pacific\nOcean", fontface="bold", size=5)
  
## Inset map
zoomrect <- data.frame(x1=c(-2130000), y1=c(-1400000), x2=c(-1700000), y2=c(-900000), r=c(1))
inset <- ggplotGrob(ggplot() + mp + statesb + outline + theme_classic() + coord_cartesian(xlim=c(-4000000,3000000), ylim=c(-3500000,2500000)) +
  geom_rect(data=zoomrect, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="black", fill=NA, size=1) +
  theme(axis.title.x=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
        axis.title.y=element_blank(), panel.border=element_rect(color="black", fill=NA, size=0.5)) +
  annotate(geom="text", x=0, y=-600000, label="North\nAmerica", size=4) +
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))) 

## Panel B
cols <- c("Sacramento Valley"="gray20","Rice Fields"="gray40", "Recovery Locations"="black")
fills <- c("Sacramento Valley"="gray80", "Rice Fields"="gray40", "Recovery Locations"="white")
recov <- coords[coords$group=="Recovery Locations",]

panelB <- ggplot() + caneor10mb + ca10mb + theme_classic() + coord_cartesian(xlim=c(-2130000,-1730000), ylim=c(-1400000,-900000)) +
  geom_polygon(data=sv, aes(x=long, y=lat, fill="Sacramento Valley", color="Sacramento Valley", group=group)) +
  geom_polygon(data=rice, aes(x=long, y=lat, fill="Rice Fields", color="Rice Fields", group=group)) +
  geom_point(data=recov, aes(x=x, y=y, shape="Recovery Locations"), fill="white", size=3) +
  theme(axis.title.x=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
        axis.title.y=element_blank(), panel.border=element_rect(color="black", fill=NA, size=0.5),
        legend.position=c(0.01, 0.3), legend.justification=c(0,1), legend.title=element_blank(),
        legend.background=element_rect(fill=NA), legend.text=element_text(size=12), legend.spacing.y=unit(-0.2, "cm")) +
  scale_shape_manual(name="", values=21) +
  scale_fill_manual(name="", values=fills) +
  scale_color_manual(name="", values=cols) +
  annotation_scale(aes(location="tr"), plot_unit="m", text_cex=1) +
  annotate(geom="text", x=-1980000, y=-1180000, label="California", size=5.5) +
  annotate(geom="text", x=-2130000, y=-890000, label="(b)", fontface="bold", size=6) +
  annotation_custom(grob=inset, xmin=-2153000, xmax=-1953000, ymin=-1065000, ymax=-920000)

gA <- ggplotGrob(panelA)
gB <- ggplotGrob(panelB)

final <- arrangeGrob(gA, gB, layout_matrix=cbind(c(1,1), c(2,2)))
grid.arrange(final)
