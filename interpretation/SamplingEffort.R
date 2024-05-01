---
title: "Sampling coverage and effort"
author: "Anna Drake"
date: "2024-04-29"
---
  
# This script examines survey coverage and effort within the projection region of the BAM species distribution models. 
# We examine coverage and effort both spatially and temporally. It is the same script used to produce an R markdown document of the same name.

### Load libraries ----------------

require(dplyr)
require(yaImpute)
require(rasterVis)
require(tidyverse)
require(terra)
require(sf)

#1. Setup ------------------------

# Colour ramp ----
col5 <- colorRampPalette(c("darkblue", "blue", "lightblue","lightgrey","beige"))

# NAD83(NSRS2007)/Conus Albers projection (epsg:5072) -----
crs<-"+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"

# Set root path for data on Google Drive ----
root<-"G:/Shared drives/BAM_NationalModels/NationalModels5.0" #PC link

# Get region shapefile ----
bcr<- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel.shp")) %>% 
  dplyr::filter(subUnit!=1) %>% # remove BCR 1
  st_transform(crs=crs) 

# Load data package ----
load(file.path(root, "Data", "04_NM5.0_data_stratify.R"))
rm("offsets","gridlist","birdlist","cov")
gc()

# Get any prediction layer to act as a blank raster ----
Canada <- readxl::read_excel(file.path(root,"NationalModels_V5_VariableList.xlsx"), sheet = "ExtractionLookup")

rast<-area<-terra::rast(file.path(root, Canada$PredictPath[1],"ERAMAP_1km_1999.tif")) %>% 
  mask(.,bcr)

rast[]<-0 # blank raster
area[]<-1 # all area of interest

# Get all sample points ----

pnts<-visit%>%.[11:12]%>% #sample lon, lat
  st_as_sf(., coords = c("lon", "lat"),crs = 4326)%>% 
  st_transform(crs=crs) %>%
  vect(.)

locations <- rasterize(pnts,rast,value=1) # put points in raster
locations[is.na(locations)]<-0 #
locations<-locations%>% mask(.,bcr) 

#plot(locations, main="Survey points: all years", legend=FALSE)


#2. Calculate the minimum distance between any survey location in the study region and any spatial point --------

tblPnts <- as.data.frame(locations, xy=TRUE) %>% subset(.,.[3]>0) #get all survey locations
colnames(tblPnts)[3] <- 'prev'
tblReg <- as.data.frame(area,xy=TRUE) #all locations
colnames(tblReg)[3] <-  'prev'

p.xy <- mutate(tblPnts, pixelID = 1:nrow(tblPnts)) %>%
  dplyr::select(pixelID, x, y, prev)
r.xy <- mutate(tblReg, pixelID = 1:nrow(tblReg)) %>% 
  dplyr::select(pixelID, x, y, prev)

p.xy2 <- p.xy %>% dplyr::select(1:3) %>% as.matrix()
r.xy2 <- r.xy %>% dplyr::select(1:3) %>% as.matrix()

if (nrow(r.xy) > 0) { d.ann <- as.data.frame(ann(
  as.matrix(p.xy2[, -1, drop = FALSE]),
  as.matrix(r.xy2[, -1, drop = FALSE]),
  k = 1,
  verbose = F)$knnIndexDist) # shortest distance in meters

d1b <- as.data.frame(cbind(r.xy2, round(sqrt(d.ann[, 2])/1000))) #change to km
names(d1b) <- c("ID", "X", "Y", "distance")

} else {
  print(spec[i])
}

r.xy <- as.data.frame(r.xy)
colnames(r.xy) <- c('ID', 'X', 'Y', 'prev')
r.xy <- as_tibble(r.xy)
d1b <- left_join(r.xy, d1b, by = c('ID', 'X', 'Y'))
final <- rast(d1b[c(2, 3, 5)], crs=crs) 
final<-final %>% crop(., bcr) %>% mask(.,bcr)

# max is 613 km - set anything above 500km to 500km
final[final>500]<-500
#writeRaster(final,"km_to_survey_point.tif", overwrite=TRUE)

plot(final, col=col5(n=10), main="Minimum distance to a survey (km)")
rm("tblPnts","tblReg","d1b")
gc()

#3. Break the data down in 5-year intervals from 1985-2020 -------------------------------
# data for each period straddles the focal year (i.e. 1985 represents data from 1983-1987). 

Out<-list()

for(i in c(1983,1988,1993,1998,2003,2008,2013,2018)){ # run each period
  # isolate 5 year period
  locSub<-subset(visit, visit$year>=i & visit$year<=i+4)
  
  locSub<-locSub%>%.[11:12]%>% #sample lon, lat
    st_as_sf(., coords = c("lon", "lat"),crs = 4326)%>% 
    st_transform(crs=crs) %>%
    vect(.)
  
  locSub <- rasterize(locSub,rast,value=1) # put points in raster with correct projection
  locSub[is.na(locSub[])]<-0 #
  locSub<-locSub%>% mask(.,bcr) 
  
  tblPnts <- as.data.frame(locSub, xy=TRUE) %>% subset(.,.[3]>0) #survey locations
  colnames(tblPnts)[3] <- 'prev'
  tblReg <- as.data.frame(area,xy=TRUE) #all locations
  colnames(tblReg)[3] <-  'prev'
  
  p.xy <- mutate(tblPnts, pixelID = 1:nrow(tblPnts)) %>%
    dplyr::select(pixelID, x, y, prev)
  r.xy <- mutate(tblReg, pixelID = 1:nrow(tblReg)) %>% 
    dplyr::select(pixelID, x, y, prev)
  
  p.xy2 <- p.xy %>% dplyr::select(1:3) %>% as.matrix()
  r.xy2 <- r.xy %>% dplyr::select(1:3) %>% as.matrix()
  
  if (nrow(r.xy) > 0) { d.ann <- as.data.frame(ann(
    as.matrix(p.xy2[, -1, drop = FALSE]),
    as.matrix(r.xy2[, -1, drop = FALSE]),
    k = 1,
    verbose = F)$knnIndexDist) # shortest distance in meters
  
  d1b <- as.data.frame(cbind(r.xy2, round(sqrt(d.ann[, 2])/1000))) #change to km
  names(d1b) <- c("ID", "X", "Y", "distance")
  
  } else {
    print(spec[i])
  }
  
  r.xy <- as.data.frame(r.xy)
  colnames(r.xy) <- c('ID', 'X', 'Y', 'prev')
  r.xy <- as_tibble(r.xy)
  d1b <- left_join(r.xy, d1b, by = c('ID', 'X', 'Y'))
  final <- rast(d1b[c(2, 3, 5)], crs=crs) 
  Out[[i]]<-final %>% crop(., bcr) %>% mask(.,bcr)
  
  #Set anything above 500km to 500km
  Out[[i]][Out[[i]][]>500]<-500
  
  rm("tblPnts","tblReg","d1b")
gc()
} #end of periods

par(mfrow=c(2,4))
for(i in c(1983,1988,1993,1998,2003,2008,2013,2018)){
plot(Out[[i]], range=c(0,500), all_levels=T, col=col5(n=10), main=paste(i+2,"km to survey", sep=" "))
}

#4. Sampling intensity: all years ---------------
# Determine sampling effort within 250 km radius of any given point

Samp<-visit%>%.[11:12]%>% 
  st_as_sf(., coords = c("lon", "lat"),crs = 4326)%>% 
  st_transform(crs=crs) %>% 
  st_buffer(.,250000)

Samp<- rasterize(Samp,rast, fun="sum") # summarize overlapping surveys
Samp[is.na(Samp)]<-(-1)
Samp<-mask(Samp,bcr)
#writeRaster(Samp,"Surveys_within_250km.tif", overwrite=TRUE)

# Bin for easier interpretation
breaks <- c(-1,0,1,5,10,50,100,500,1000,5000,10000,50000,100000)

# Create raster layers with values binned by powers of 10
rc2 <- classify(Samp, breaks, include.lowest=TRUE, brackets=TRUE)
rm("Samp")
gc()

# Plot the binned raster
plot(rc2, col=rev(col5(n=12)), main="Survey intensity within 250km radius: all data")

#5. Break data into 5 year intervals from 1985-2020 to get changes in temporal coverage ---------------
# Contributing data straddles the focal years, as above (1985= 1983-1987). 
# Assuming habitat/climate/disturbance relationships are stable over time, temporal 
# coverage matters most for the reliability of "trend" ("year") in the model. "-1" == no coverage 

Out<-list()

for (i in c(1983,1988,1993,1998,2003,2008,2013,2018)){ # run each period
  # isolate 5 year period
  SubSamp<-subset(visit, visit$year>=i & visit$year<=i+4 )  
  
  SubSamp<-SubSamp%>%.[11:12]%>% 
    st_as_sf(., coords = c("lon", "lat"),crs = 4326)%>% 
    st_transform(crs=crs) %>% 
    st_buffer(.,250000)
  
  #sum in raster  
  SubSamp<- rasterize(SubSamp,rast, fun="sum") # summarize overlapping surveys
  SubSamp[is.na(SubSamp[])]<-(-1)
  SubSamp<-mask(SubSamp,bcr)
  
  #bin values
  rc2 <- classify(SubSamp, breaks, include.lowest=TRUE, brackets=TRUE)
  Out[[i]]<-rc2
  rm("SubSamp")
gc()
} #end of periods

par(mfrow=c(2,4))
for(i in c(1983,1988,1993,1998,2003,2008,2013,2018)){
 plot(Out[[i]], range=c(-1,100000), all_levels=T, col=rev(col5(n=13)), main=paste(i+2,"coverage", sep=" "), legend=F)
  }
  
################## END OF CODE ############################
