library(raster)
library(sf)
library(maptools)
library(dplyr)
library(data.table)
library(reshape2)
library(raster)
library(geosphere)
library(sp)
library(rgdal)
require(rgeos)
library(MODISTools)
library(tidyr)
library(RODBC)
library(geojsonio)
library(rgee)
ee_install_set_pyenv(
  py_path = "C:/Users/mehou10/AppData/Local/r-miniconda", # Change it for your own Python PATH
  py_env = "rgee" # Change it for your own Python ENV
)
ee_check()
ee_Initialize()

#Set params
w <- "./data"
out <- file.path(w,"outsample")
wt_bind <- "./01_wtBind.RData"

## Load WT data
load(wt_bind)


# Proj
LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
LAEA <- CRS("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
DD <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")  #EPSG:4386
TM10 <- CRS("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +datum=NAD83 +units=m +no_defs")
LCC49 <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=49 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")


# Reproject location 
#SSBBS_DD <- SSBBS
#coordinates(SSBBS_DD) <- c("longitude", "latitude")
#proj4string(SSBBS_DD) <- DD



#---------------------------------------------------------------
#######    TEMP STEP: Work with sample n= 1000
#---------------------------------------------------------------
location <-unique(wt_bind[c("organization", "project", "sensor", "status", "location", "latitude", "longitude")])
PClocation <-PClocation[sample(nrow(PClocation), 10000), ]
PKEYcombo <-subset(PKEYcombo, PKEYcombo$SS %in% PClocation$SS)


#---------------------------------------------------------------

PClocation_DD <- PClocation
PClocation_DD <- na.omit(PClocation_DD)
coordinates(PClocation_DD) <- c("longitude", "latitude")
proj4string(PClocation_DD) <- DD

PClocation_LCC <- spTransform(PClocation_DD,LCC)
PClocationLCC <- as.data.frame(PClocation_LCC)
PClocation_LAEA <- spTransform(PClocation_DD,LAEA)
PClocationLAEA <- as.data.frame(PClocation_LAEA)
PClocation_TM10 <- spTransform(PClocation_DD,TM10)
PClocationTM10 <- as.data.frame(PClocation_TM10)

PClocation_LCC49 <- spTransform(PClocation_DD,LCC49)
PClocationLCC49 <- as.data.frame(PClocation_LCC49)

#Reproject visit using location
#PKEYcombo_DD <- left_join(PKEYcombo, PClocation, by="SS")
PKEYcombo_DD <- merge(PKEYcombo, PClocation, by="SS")

PKEYcombo_DD <- na.omit(PKEYcombo_DD)
PKEYcomboDD <- as.data.frame(PKEYcombo_DD)
coordinates(PKEYcombo_DD) <- c("longitude", "latitude")
proj4string(PKEYcombo_DD) <- DD

PKEYcombo_TM10 <- left_join(PKEYcombo, PClocationTM10, by="SS")
PKEYcombo_TM10 <- na.omit(PKEYcombo_TM10)
PKEYcomboTM10 <- as.data.frame(PKEYcombo_TM10)
coordinates(PKEYcombo_TM10) <- c("longitude", "latitude")
proj4string(PKEYcombo_TM10) <- TM10

PKEYcombo_LCC <- left_join(PKEYcombo, PClocationLCC, by="SS")
PKEYcombo_LCC <- na.omit(PKEYcombo_LCC)
PKEYcomboLCC <- as.data.frame(PKEYcombo_LCC)
coordinates(PKEYcombo_LCC) <- c("longitude", "latitude")
proj4string(PKEYcombo_LCC) <- LCC

PKEYcombo_LAEA <- left_join(PKEYcombo, PClocationLAEA, by="SS")
PKEYcombo_LAEA <- na.omit(PKEYcombo_LAEA)
PKEYcomboLAEA <- as.data.frame(PKEYcombo_LAEA)
coordinates(PKEYcombo_LAEA) <- c("longitude", "latitude")
proj4string(PKEYcombo_LAEA) <- LAEA

PKEYcombo_LCC49 <- left_join(PKEYcombo, PClocationLCC49, by="SS")
PKEYcombo_LCC49 <- na.omit(PKEYcombo_LCC49)
PKEYcomboLCC49 <- as.data.frame(PKEYcombo_LCC49)
coordinates(PKEYcombo_LCC49) <- c("longitude", "latitude")
proj4string(PKEYcombo_LCC49) <- LCC49
#---------------------------------------------------------------
 #          TEMP    CONSTRAINT
#---------------------------------------------------------------
###### CHECK if NA and blank
apply(PKEYcombo, 2, function(x) any(is.na(x))) #NAs not present
apply(PKEYcombo, 2, function(x) any(x=="")) #blanks not present

apply(PClocation, 2, function(x) any(is.na(x))) #NAs present lat long column
apply(PClocation, 2, function(x) any(x==""))  

PClocation <-PClocation[!(is.na(PClocation$latitude)),]
#dataBLANK <- PClocation[PClocation$latitude=="",]

#check if lat long are NA in db (not in location tbl. DON'T KNOW WHY!!!!!)
PClocation2 <- sqlFetch(con, "location")
PClocation2[PClocation2$SS_V4=="BCCA:08MM10:306139",]
nrow(PClocation2) #149912
PClocation2 <- PClocation2[,c(3,5:6)]
names(PClocation2) <- c("SS", "latitude", "longitude")
PClocation3 <-  unique(PClocation2[ , ])
nrow(PClocation3) #149912 good all unique.
#Now check if there are PKEY without any location and vice versa

PKEYcombo2 <- sqlFetch(con, "pc_visit")
PKEYlocation <- unique(PKEYcombo2$location_name_4)
length(PKEYlocation) #146328
subset(PKEYlocation, !(PKEYlocation %in% PClocation3$SS)) # have pkey but no location "MAVEXP:A:a28" "MAVEXP:A:A30"
subset(PClocation3, !(PClocation3$SS %in% PKEYlocation)) # 21500 have location but no pkey
PKEYcombo2[PKEYcombo2$location_name_4=="BMS:TLU-074:36",]
#---------------------------------------------------------------
#---------------------------------------------------------------




####### Fixed layers ########
#Landform (100-m)
#-------------------------------------------------
#  
#                 stack Topo
#-------------------------------------------------
TPI <- raster(file.path(w,"topoedaphic/TPI.tif"))
TRI <- raster(file.path(w,"topoedaphic/TRI.tif"))
slope <- raster(file.path(w,"topoedaphic/slope.tif"))
roughness <- raster(file.path(w,"topoedaphic/roughness.tif"))
topo <- stack(TPI,TRI,slope,roughness)
PClocation_LAEA <- cbind(PClocation_LAEA,raster::extract(topo,PClocation_LAEA))
#-------------------------------------------------

#-------------------------------------------------
#               Road
usroaddist <- raster(file.path(w,"transportation/USroaddistLAZEA.tif"))
canroaddist <- raster(file.path(w,"transportation/canroaddist1.tif"))
PClocation_LAEA <- cbind(PClocation_LAEA,raster::extract(canroaddist,PClocation_LAEA))
names(PClocation_LAEA)[ncol(PClocation_LAEA)] <- "CanRoadDist"
PClocation_LAEA <- cbind(PClocation_LAEA,raster::extract(usroaddist,PClocation_LAEA))
names(PClocation_LAEA)[ncol(PClocation_LAEA)] <- "USRoadDist"
countries <- shapefile(file.path(w,"basemaps/NorthAmericaLazea.shp"))
countries <- spTransform(countries,LAEA)
PClocation_LAEA <- cbind(PClocation_LAEA,raster::extract(countries,PClocation_LAEA))
PClocation_LAEA <- PClocation_LAEA[,c(1:7,10)]
names(PClocation_LAEA)[8] <- "Country"
PClocation_LAEA$RoadDist <- pmin(PClocation_LAEA$CanRoadDist,PClocation_LAEA$USRoadDist, na.rm=TRUE)

write.csv(PClocation_LAEA,file=file.path(out,"SS_TopoRoad.csv"),row.names=FALSE)

#-------------------------------------------------
lf <- raster(file.path(w,"topoedaphic/lf_lcc1.tif"))
PClocation_LCC <- cbind(PClocation_LCC,raster::extract(lf,PClocation_LCC))
names(PClocation_LCC)[ncol(PClocation_LCC)] <- "lf"

#Adaptwest baseline climate variables
cur <- file.path(w,"climate/baseline19812010/")
setwd(cur)
clim <- list.files(cur, pattern =".asc$")
curclim<-stack(clim)
PClocation_LCC <- cbind(PClocation_LCC,raster::extract(curclim,PClocation_LCC))

write.csv(PClocation_LCC,file=file.path(out,"SS_LandformClimate.csv"),row.names=FALSE)

#YearTier1 <- c(2001:2011)
#YearTier2 <- c(1991:2000,2012:2020)
#-------------------------------------------------

####### Biomass layers ######
#  
#-------------------------------------------------
###kNN biomass layers, 2001 (250-m)
  #b2001 <- list.files(file.path(w,"biomass/2001"),pattern="tif$")
  #setwd(file.path(w,"biomass/2001"))
  #bs2001 <- stack(raster(b2001[1]))
  #for (i in 2:length(b2001)) {bs2001 <- addLayer(bs2001, raster(b2001[i]))}
  #names(bs2001) <- gsub("NFI_MODIS250m_2001_kNN_","",names(bs2001))
  #writeRaster(bs2001, file=file.path(w,"biomass/bs2001_250"))
bs2001 <- brick(file.path(w,"biomass/bs2001_250.grd" ))
PClocation_LCC2001 <- cbind(PClocation_LCC,raster::extract(bs2001,PClocation_LCC))

#bs2001_Gauss750 <- brick("G:/Boreal/NationalModelsV2/bs2001_750.grd")
bs2001_Gauss750 <- brick(file.path(w, "biomass/bs2001_750.grd"))
names(bs2001_Gauss750) <- gsub("SpeciesGroups","Landsc750",names(bs2001_Gauss750))
names(bs2001_Gauss750) <- gsub("Species","Landsc750",names(bs2001_Gauss750))
names(bs2001_Gauss750) <- gsub("Structure","Landsc750",names(bs2001_Gauss750))
names(bs2001_Gauss750) <- gsub("LandCover","Landsc750",names(bs2001_Gauss750))
PClocation_LCC2001 <- cbind(PClocation_LCC2001,raster::extract(bs2001_Gauss750,PClocation_LCC2001))



###kNN biomass layers, 2011 (250-m)
  #b2011 <- list.files("E:/GIS/landcover/Beaudoin/Processed_sppBiomass/2011/",pattern="tif$")
  #setwd("E:/GIS/landcover/Beaudoin/Processed_sppBiomass/2011/")
  #bs2011 <- stack(raster(b2011[1]))
  #for (i in 2:length(b2011)) {bs2011 <- addLayer(bs2011, raster(b2011[i]))}
  #names(bs2011) <- gsub("NFI_MODIS250m_2011_kNN_","",names(bs2011))
bs2011 <- brick(file.path(w,"biomass/bs2011_250m.grd" ))
PClocation_LCC2011 <- cbind(PClocation_LCC,raster::extract(bs2011,PClocation_LCC))

bs2011_Gauss750 <- brick(file.path(w,"biomass/bs2011_750.grd"))
names(bs2011_Gauss750) <- gsub("SpeciesGroups","Landsc750",names(bs2011_Gauss750))
names(bs2011_Gauss750) <- gsub("Species","Landsc750",names(bs2011_Gauss750))
names(bs2011_Gauss750) <- gsub("Structure","Landsc750",names(bs2011_Gauss750))
names(bs2011_Gauss750) <- gsub("LandCover","Landsc750",names(bs2011_Gauss750))
PClocation_LCC2011 <- cbind(PClocation_LCC2011,raster::extract(bs2011_Gauss750,PClocation_LCC2011))

PKEYcombo_2001 <- left_join(PKEYcombo[PKEYcombo$YEAR<2006,], as.data.frame(PClocation_LCC2001), by="SS")
PKEYcombo_2011 <- left_join(PKEYcombo[PKEYcombo$YEAR>2005,], as.data.frame(PClocation_LCC2011), by="SS")
PKEYcombo_dat <- rbind(PKEYcombo_2001,PKEYcombo_2011)
write.csv(PKEYcombo_dat,file=file.path(out,"PKEYcombo_VegBiomass.csv"),row.names=FALSE)


####### Annual layers ######

#McKenney annual climate variables
annclim <- function(surveys, yr) {
  cmi <- raster(paste0(w,"/climate/NRCAN/Annual/",yr,"/cmi60_sum.asc"))
  clim <- stack(lapply(grep("bio",list.files(paste0(w,"/climate/NRCAN/Annual/",yr),pattern=".tif$",full.names=TRUE),value=TRUE),FUN=raster))
  cmi <- resample(cmi,clim[[1]])
  clim <- addLayer(clim,cmi)
  data <- cbind(surveys,raster::extract(clim,surveys))
  return(data)
}

CanadaMODISLC <- function (surveys, yr) {
  nalcx <- raster(paste0(w,"/landcover/LCTS_2000-2011/LCTS19c",yr,".tif"))
  data <- cbind(surveys,raster::extract(nalcx,surveys))
  names(data)[ncol(data)]<-"LCTS250m"
  return(data)
}

MODIS500 <- function(surveys,yr){
  modx <- raster(paste0(w,"/landcover/LCCMODIS500/LC500_",yr,".tif"))  
  data <- cbind(surveys,raster::extract(modx,surveys))
  names(data)[ncol(data)]<-"LC500m"
  return(data)
}

# MODIS500 <- function (surveys, yr) {
#   df <- surveys
#   df$site_name = df$SS
#   df <- df[,c(6,5,4)]
#   names(df) <- c("site_name","lon","lat")
#   lcdata <- mt_batch_subset(df = df,
#                             product = "MCD12Q1",
#                             band = "LC_Type1",
#                             start = paste0(yr,"-01-01"),
#                             end = paste0(yr,"-12-31"),
#                             out_dir = "E:/GIS/landcover/LCCMODIS500",
#                             internal = TRUE)
#   names(lcdata)[ncol(lcdata)]<-"LCTS500m"
#   return(lcdata)
# }

b1 <- brick(file.path(w,"landcover/NASA_ABoVE/1_to_17_mosaic.tif"))
ABOVE <- function(surveys, yr) {
  data <- cbind(surveys,raster::extract(b1[[yr-1983]],surveys))
  names(data)[ncol(data)] <- "ABOVELC"
  return(data)
}

VLCE25 <- function(surveys, yr) {
  hermo <- ee$ImageCollection("projects/sat-io/open-datasets/CA_FOREST_LC_VLCE2")$ # read the data
    select('b1')$ # select landcover
    filter(ee$Filter$calendarRange(yr, yr, "year"))$
    toBands()# filter the year
  temp <- ee_extract(
    x = hermo,
    y = surveys["PKEY"],
    scale = 25,
    sf = FALSE
  )
  if(ncol(temp)<2){
    temp[,2] <- NA
  }
  yr1dat_df <- cbind(surveys,temp[,2])
  yr1dat_df <- yr1dat_df %>% st_drop_geometry()
  colnames(yr1dat_df)[4] <- "VLCE25"
  return(yr1dat_df)
}
# Extract NLCD (USA)
NLCD_con <- function(surveys, yr) {
  if (yr<2002) { 
    yr <- 2001
  } else if (2001 < yr & yr < 2005) {
    yr <- 2004
  } else if (2004 < yr & yr < 2009) {
    yr <- 2008
  } else if (2008 < yr & yr < 2012) {
    yr <- 2011
  } else if (2011 < yr & yr < 2014) {
    yr <- 2013
  } else if (2013 < yr & yr < 2017) {
    yr <- 2016
  } else {
    yr <- 2019
  }  
  nlcd <- ee$ImageCollection("USGS/NLCD_RELEASES/2019_REL/NLCD")$ # read the data
    select('landcover')$ # select landcover
    filter(ee$Filter$eq('system:index', as.character(yr)))
  temp <- ee_extract(
    x = nlcd,
    y = surveys["PKEY"],
    scale = 30,
    sf = FALSE
  )
  if(ncol(temp)<2){
    temp[,2] <- NA
  }
  yr1dat_df <- cbind(surveys,temp[,2])
  #yr1dat_df <- yr1dat_df %>% st_drop_geometry()
  colnames(yr1dat_df)[5] <- "NLCD_con"
  return(yr1dat_df)
}
NLCD_ak <- function(surveys, yr) {
  if (yr<2006) { 
    yr <- 2001
  } else if (2005 < yr & yr < 2014) {
    yr <- 2011
  } else {
    yr <- 2016
  }  
  nlcd <- ee$ImageCollection("USGS/NLCD_RELEASES/2016_REL")$ # read the data
    select('landcover')$ # select landcover
    filter(ee$Filter$eq('system:index', paste0(as.character(yr), "_AK")))
  temp <- ee_extract(
    x = nlcd,
    y = surveys["PKEY"],
    scale = 30,
    sf = FALSE
  )
  if(ncol(temp)<2){
    temp[,2] <- NA
  }
  yr1dat_df <- cbind(surveys,temp[,2])
  #yr1dat_df <- yr1dat_df %>% st_drop_geometry()
  colnames(yr1dat_df)[6] <- "NLCD_ak"
  
  return(yr1dat_df)
}

# Extract Global Canopy Height 10m
GCH20 <- function(surveys) {
  nlang <- ee$Image("users/nlang/ETH_GlobalCanopyHeight_2020_10m_v1")$ # read the data
    select('b1') 
  temp <- ee_extract(
    x = nlang,
    y = surveys["PKEY"],
    scale = 25,
    sf = FALSE
  )
  if(ncol(temp)<2){
    temp[,2] <- NA
  }
  yr1dat_df <- cbind(surveys,temp[,2])
  yr1dat_df <- yr1dat_df %>% st_drop_geometry()
  colnames(yr1dat_df)[4] <- "GCH20_height"
  return(yr1dat_df)
}
###############################

yearlist <- sort(unique(PKEYcombo$YEAR))

yearlist1 <- yearlist[yearlist < 2015]
yearlist1 <- yearlist1[yearlist1 > 1983]
yrdat <- setNames(data.frame(matrix(ncol = 22, nrow = 0)), c(names(PKEYcombo_TM10),"ABOVELC")) 
for (i in 1:length(yearlist1)) {
  yr1 <- PKEYcombo_TM10[PKEYcombo_TM10$YEAR==yearlist1[i],]
  yr1dat <- ABOVE(yr1, as.integer(yearlist1[i]))
  yrdat <- rbind(yrdat,yr1dat@data)
}
write.csv(yrdat,file=paste0(out,"/PKEYannualLCABOVE.csv"),row.names=FALSE)

clim <- stack(lapply(grep("bio",list.files(paste0(w, "/climate/NRCAN/Annual/",2005),pattern=".tif$",full.names=TRUE),value=TRUE),FUN=raster))
cmi <- raster(paste0(w, "/climate/NRCAN/Annual/2005/cmi60_sum.asc"))
cmi <- resample(cmi,clim[[1]])
clim <- addLayer(clim,cmi)
yearlist2 <- yearlist[yearlist < 2019]
yrdatclim <- setNames(data.frame(matrix(ncol = 41, nrow = 0)), c(names(PKEYcombo_DD),names(clim))) 
for (i in 1:length(yearlist2)) {
  yr1 <- PKEYcombo_DD[PKEYcombo_DD$YEAR==yearlist2[i],]
  yr1dat <- annclim(yr1, as.integer(yearlist2[i]))
  yrdatclim <- rbind(yrdatclim,as.data.frame(yr1dat))
}
write.csv(yrdatclim,file=paste0(out,"/PKEYannualclim.csv"),row.names=FALSE)

yearlist3 <- yearlist[yearlist < 2012]
yearlist3 <- yearlist3[yearlist3 > 1999]
yrdatlc <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c(names(PKEYcombo_LCC),"LCTS250m")) 
for (i in 1:length(yearlist3)) {
  yr1 <- PKEYcombo_LCC[PKEYcombo_LCC$YEAR==yearlist3[i],]
  yr1dat <- CanadaMODISLC(yr1, as.integer(yearlist3[i]))
  yrdatlc <- rbind(yrdatlc,yr1dat@data)
}
write.csv(yrdatlc,file=paste0(out,"/PKEYannualLC250.csv"),row.names=FALSE)

yearlist4 <- yearlist[yearlist < 2020]
yearlist4 <- yearlist4[yearlist4 > 2000]
yrdatlc <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c(names(PKEYcombo_DD),"LC500m")) 
for (i in 1:length(yearlist4)) {
  yr1 <- PKEYcombo_DD[PKEYcombo_DD$YEAR==yearlist4[i],]
  yr1dat <- MODIS500(yr1, as.integer(yearlist4[i]))
  yrdatlc <- rbind(yrdatlc,yr1dat@data)
}
write.csv(yrdatlc,file=paste0(out,"/PKEYannualLC500.csv"),row.names=FALSE)

# VLCE25
yearlist <- sort(unique(PKEYcombo_DD$YEAR))
yearlist1 <- yearlist[yearlist < 2020]
yearlist1 <- yearlist1[yearlist1 > 1983]
yrdat <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c(names(PKEYcombo_DD),"VLCE25")) 
for (i in 1:length(yearlist1)) {
  m <- PKEYcombo_DD[PKEYcombo_DD$YEAR==yearlist1[i],]
  n <- nrow(m)
  r  <- rep(1:ceiling(n / 5000), each = 5000, length.out = nrow(m) )[1:n]
  d <- split(m, r) # split df to list
  
  for(j in 1:length(d)) {
    PKEYcombo_DD_sf <- st_as_sf(d[[j]])
    PKEYextract <- VLCE25(PKEYcombo_DD_sf, as.integer(yearlist1[i]))
    yrdat <- rbind(yrdat, PKEYextract[,1:4])
  }
}
write.csv(yrdat,file=paste0(out,"/PKEYannualVLCE25.csv"),row.names=FALSE)

# NLCD
  # Process only US to optimize extraction
countries <- shapefile(file.path(w,"basemaps/NorthAmericaLazea_v2.shp"))
countries <- spTransform(countries,DD)
PKEYcombo_DD_1 <- cbind(PKEYcombo_DD,raster::extract(countries,PKEYcombo_DD))
PKEYcombo_DD_1 <- PKEYcombo_DD_1[,c(1:3,8)]
names(PKEYcombo_DD_1)[4] <- "Country"

PKEYcombo_DD_US <- subset(PKEYcombo_DD_1, Country!="CAN" | is.na(Country))

yearlist <- sort(unique(PKEYcombo_DD_US$YEAR))
yearlist1 <- yearlist[yearlist < 2020]
yearlist1 <- yearlist1[yearlist1 > 2000]
yrdat <- setNames(data.frame(matrix(ncol =8, nrow = 0)), c(names(PKEYcombo_DD_US), "longitude", "latitude", "NLCD_con", "NLCD_ak"))
for (i in 1:length(yearlist1)) {
  m <- PKEYcombo_DD_US[PKEYcombo_DD_US$YEAR==yearlist1[i],]
  PKEYcombo_DD_sf <- st_as_sf(m)
  PKEYextract_con <- nlcd_con(PKEYcombo_DD_sf, as.integer(yearlist1[i]))
  print(nrow(PKEYextract_con))
  PKEYextract_ak <- nlcd_ak(PKEYextract_con, as.integer(yearlist1[i]))
  print(nrow(PKEYextract_ak))
  PKEYextract <-as.data.frame(as(PKEYextract_ak, Class="Spatial"))
  colnames(PKEYextract)[7] <- "longitude"
  colnames(PKEYextract)[8] <- "latitude"
  yrdat <- rbind.fill(yrdat, PKEYextract)
}

write.csv(yrdat,file=paste0(out,"/PKEYannualNLCD.csv"),row.names=FALSE)

# GCH2020
yearlist <- sort(unique(PKEYcombo_DD$YEAR))
yearlist1 <- yearlist[yearlist > 2000]
yrdat <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c(names(PKEYcombo_DD),"GCH20_height")) 
for (i in 1:length(yearlist1)) {
  m <- PKEYcombo_DD[PKEYcombo_DD$YEAR==yearlist1[i],]
  n <- nrow(m)
  r  <- rep(1:ceiling(n / 5000), each = 5000, length.out = nrow(m) )[1:n]
  d <- split(m, r) # split df to list
  
  for(j in 1:length(d)) {
    PKEYcombo_DD_sf <- st_as_sf(d[[j]])
    PKEYextract <- GCH20(PKEYcombo_DD_sf)
    yr1dat_df <- yr1dat %>% st_drop_geometry()
    yrdat <- rbind(yrdat, PKEYextract[,1:4])
  }
}

write.csv(yrdat,file=paste0(out,"/PKEYannualGCH20.csv"),row.names=FALSE)

# LANDFIRE
LF <- file.path(w,"landcover/LandFire")
  ##Stack raster
  #LFak <- list.files(file.path(LF,"AK"),pattern="tif$")
  #setwd(file.path(LF,"AK"))
  #temp<-raster(LFak[1])
  #LF_ak <- stack(temp)
  #for (i in 2:length(LFak)) {
  #  LF_ak <- addLayer(LF_ak, raster(LFak[i]))}
  #writeRaster(LF_ak, file=file.path(LF,"LF_AK_30"))

  #LFcon <- list.files(file.path(LF,"CONUS"),pattern="tif$")
  #setwd(file.path(LF,"CONUS"))
  #temp<-raster(LFcon[1])
  #LF_con <- stack(temp)
  #for (i in 2:length(LFcon)) {
  #  LF_ak <- addLayer(LF_con, raster(LFcon[i]))}

  #writeRaster(LF_con, file=file.path(LF,"LF_CON_30"))


### Combine Data ####
toproad <- read.csv(file.path(out,"SS_TopoRoad.csv"))
lfclimate <- read.csv(file.path(out,"SS_LandformClimate.csv"))
vegbiomass <- read.csv(file.path(out,"PKEYcombo_VegBiomass.csv"))
annualclim <- read.csv(file.path(out,"PKEYannualclim.csv"))
annualABOVE <- read.csv(file.path(out,"PKEYannualLCABOVE.csv"))
annualMODIS250 <- read.csv(file.path(out,"PKEYannualLC250.csv"))
annualMODIS500 <- read.csv(file.path(out,"PKEYannualLC500.csv"))
annualVLCE25 <- read.csv(file.path(out,"PKEYannualVLCE25.csv"))

datcombo <- left_join(toproad[,c(1,10,11,2:9)],lfclimate[,c(1:28)],by="SS")
datcombo1 <- left_join(annualclim[,c(1:23)],vegbiomass[,c(3,31:216)],by="PKEY")
datcombo1 <- left_join(datcombo1,annualABOVE[,c(3,4)],by="PKEY")
datcombo1 <- left_join(datcombo1,annualMODIS250[,c(3,4)],by="PKEY")
datcombo1 <- left_join(datcombo1,annualMODIS500[,c(3,4)],by="PKEY")
datcombo1 <- left_join(datcombo1,annualVLCE25[,c(3,4)],by="PKEY")
datcombo2 <- right_join(datcombo,datcombo1,by="SS")
#datcombo2001 <- datcombo2[datcombo2$YEAR>2000,]

IGBP <- read.csv(file.path(w, "XWalk/IGBP_XWalk.csv"))
NALCMS <- read.csv(file.path(w, "XWalk/NALCMS_XWalk.csv"))
ABOVE <- read.csv(file.path(w, "XWalk/ABOVE_XWalk.csv"))
VLCE25 <- read.csv(file.path(w, "XWalk/VLCE_XWalk.csv"))
datcombo2 <- left_join(datcombo2,ABOVE,by = c("ABOVELC"="ABOVEClass"))
datcombo2 <- left_join(datcombo2,NALCMS,by = c("LCTS250m"="NALCMSClass"))
datcombo2 <- left_join(datcombo2,IGBP,by = c("LC500m"="IGBPClass"))
datcombo2 <- left_join(datcombo2,VLCE25,by = c("VLCE25"="VLCEClass"))
datcombo2$LCX <- ifelse(datcombo2$ABOVEX>0, datcombo2$ABOVEX, ifelse(datcombo2$VLCEX>0, datcombo2$VLCEX, ifelse(datcombo2$NALCMSX>0, datcombo2$NALCMSX, datcombo2$IGBPX)))
datcombo2$LCX <- ifelse(is.na(datcombo2$LCX), datcombo2$ABOVEX, datcombo2$LCX)
datcombo2$LCX <- ifelse(is.na(datcombo2$LCX), datcombo2$VLCEX, datcombo2$LCX)
datcombo2$LCX <- ifelse(is.na(datcombo2$LCX), datcombo2$NALCMSX, datcombo2$LCX)
datcombo2$LCX <- ifelse(is.na(datcombo2$LCX), datcombo2$IGBPX, datcombo2$LCX)
write.csv(datcombo2, file=file.path(out,"PKEYannualcombo.csv"),row.names=FALSE)

table(as.factor(datcombo2$YEAR), as.factor(datcombo2$LCX))
table(as.factor(datcombo2$YEAR), as.factor(datcombo2$LC500m))
table(as.factor(datcombo2$YEAR), as.factor(datcombo2$IGBPX))
table(as.factor(datcombo2$YEAR), as.factor(datcombo2$NALCMSX))
