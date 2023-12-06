#################################################
#  Producing Prediction Rasters for BAM NAM 4.1
#  February, 2023
#################################################

#install.packages('terra', repos='https://rspatial.r-universe.dev')
library(terra)
library(raster)

getwd()
setwd("/Users/annadrake/Downloads/")
e<-extent(-4100000, 3200000, 1673000, 6200000 )
FormatRaster<-rast("image_export.tif")  #any download from GEE for files that need more alignment 

#var region = ee.Geometry.BBox(-169.57, 71.22,-50.39, 38.09);
#Export.image.toDrive({
#  image: image,
#  description: 'image_export',
#  region: region,
#  scale: 1000,
#  maxPixels: 1e13,
#  crs: 'EPSG:5072'
#});

Format<-crop(FormatRaster,e)
values(Format)<-0

## CLIMATE NORMALS
setwd("/Users/annadrake/Downloads/ClimateNA_V7_3_1991_2020_Norm/")

##Crop and reproject and align Climate Normals
files<-list.files(pattern="Norm",all.files=TRUE, full.names=FALSE)

Norm<-list()
for (j in 1:length(files)) {
  Norm[[j]]<-rast(files[[j]])
  Norm[[j]]<-terra::project(Norm[[j]],Format,mask=TRUE,align=TRUE,method="near")  
  Norm[[j]]<-crop(Norm[[j]],e) 
writeRaster(Norm[[j]],paste("AlbersProj/",gsub("Normal_1991_2020_","",files[j]),sep=""))
}


## WETLAND VARIABLES   
WOccur<-rast("Wet_Occur.tif")
WOccur5k<-focal(WOccur,w=matrix(1,5,5),fun=mean,na.rm=TRUE)#get 5k average at 1k resolution

WRecur<-rast("Wetland_recur.tif")
WSeason<-rast("Wetland_seasonality.tif")

WOccur<-crop(WOccur, e)
WOccur5k<-crop(WOccur5k, e)
WRecur<-crop(WRecur, e)
WSeason<-crop(WSeason, e)

writeRaster(WOccur,"Wet_Occurrence.tif")
writeRaster(WOccur5k,"Wet_Occurrence_5k.tif")
writeRaster(WRecur,"Wet_Recurrence.tif")
writeRaster(WSeason,"Wet_Seasonality.tif")

## TOPOGRAPHY ####
files<-c("mTPI_Canada.tif","hli3cl.tif","TRI.tif","roughness.tif")
Topo<-list()
#for (j in 1:length(files)){
  for (j in 3:4){
  Topo[[j]]<-rast(files[[j]])
  Topo[[j]]<-terra::project(Topo[[j]],Format,mask=TRUE,align=TRUE,method="near")  
  Topo[[j]]<-crop(Topo[[j]],e) 
writeRaster(Topo[[j]],paste("Topography/",files[j],sep=""))
}

##### EVI GREENUP ######
# Days since Jan 1, 1970 - 365*x years since, to get annual value

setwd("/Users/annadrake/Downloads/EVI_Greenup/")

##Crop and reproject and align Climate Normals
Greenfiles<-list.files(pattern=".tif",all.files=TRUE, full.names=FALSE)
Period<-c(34,39,44,49,34,39,44,49)#confirm order


Green<-list()
for (j in 1:length(Greenfiles)){
  Green[[j]]<-rast(Greenfiles[[j]])
  
  Green[[j]][Green[[j]] < quantile(values(Green[[j]]),0.01, na.rm=TRUE)] <- NA  # remove the most extreme values - currently annual range is more than 500 days, so some kind of error in source files here
  Green[[j]][Green[[j]] > quantile(values(Green[[j]]),0.99, na.rm=TRUE)] <- NA
  
  #Green[[j]]<-terra::project(Green[[j]],Format,mask=TRUE,align=TRUE,method="near")  
  Green[[j]]<-crop(Green[[j]],e) 
  Green[[j]] <- Green[[j]] - Period[j]*365
  writeRaster(Green[[j]],paste("Cleaned_500m",Greenfiles[j],sep=""), overwrite=TRUE)}

summary(values(Green[[8]]))
18234-17716 #518 days in the year???
#outer quartiles fall within 36 days of each other


### Greenup average 2001-2005 #####
# 31-35 years
setwd("/Users/annadrake/Downloads/Av Greenup/")

##Crop and reproject and align Climate Normals
Greenfiles<-list.files(pattern=".tif",all.files=TRUE, full.names=FALSE)
Period<-c(31,32,33,34,35,31,32,33,34,35)#confirm order


Green<-list()
for (j in 1:length(Greenfiles)){
  Green[[j]]<-rast(Greenfiles[[j]])
  
  Green[[j]][Green[[j]] < quantile(values(Green[[j]]),0.01, na.rm=TRUE)] <- NA  # remove the most extreme values - currently annual range is more than 500 days, so some kind of error here
  Green[[j]][Green[[j]] > quantile(values(Green[[j]]),0.99, na.rm=TRUE)] <- NA
  
  Green[[j]]<-terra::project(Green[[j]],Format,mask=TRUE,align=TRUE,method="near")  
  Green[[j]]<-crop(Green[[j]],e) 
  Green[[j]] <- Green[[j]] - Period[j]*365
}
  AvDormancy<-rast(Green[c(1:5)]) %>% terra::mean(.)
  AvGreenup<-rast(Green[c(6:10)]) %>% terra::mean(.)

  #writeRaster(AvDormancy,"Cleaned_500m_Dormancy2001_2005.tif", overwrite=TRUE)
  #writeRaster(AvGreenup,"Cleaned_500m_Greenup2001_2005.tif", overwrite=TRUE)
  
  writeRaster(AvDormancy,"StandardDormancy2001_2005.tif", overwrite=TRUE)
  writeRaster(AvGreenup,"StandardGreenup2001_2005.tif", overwrite=TRUE)
  
?mean

#### HUMAN FOOTPRINT CANADA #######
HF<-rast("cum_threat2020.02.18.tif")
HF_algn<-terra::project(HF,Format,mask=TRUE,align=TRUE,method="near")  
HF_algn<-crop(HF_algn,e) 
writeRaster(HF_algn,"CAN_HF_1km.tif")

HF_algn5k<-focal(HF_algn,w=matrix(1,5,5),fun=mean,na.rm=TRUE)#get 5k average at 1k resolution
writeRaster(HF_algn5k,"CAN_HF_5km.tif")
plot(HF_algn5k)

#### Global HUMAN FOOTPRINT  ####
GHF<-rast("Global_HF_Conus.tif")
GHF<-crop(GHF,e)
GHF<-subst(GHF, 128, NA)
writeRaster(GHF,"Global_HF_Conus.tif")

GHF_algn<-terra::project(GHF,Format,mask=TRUE,align=TRUE,method="near")  
writeRaster(GHF_algn,"Global_HF_1km.tif")
G_HF_algn5k<-focal(GHF_algn,w=matrix(1,5,5),fun=mean,na.rm=TRUE)#get 5k average at 1k resolution
writeRaster(G_HF_algn5k,"Global_HF_5km.tif")
plot(G_HF_algn5k)

### PEATLAND
Peat<-rast("Peat-ML_global_peatland_extent.tif")
plot(Peat_algn)
crop<-extent(-170, -0, 0, 80)
Peat<-crop(Peat,crop)
Peat_algn<-terra::project(Peat,Format,mask=TRUE,align=TRUE,method="near")  
Peat_algn<-crop(Peat_algn,e)
writeRaster(Peat_algn,"Peatland_1km.tif")


### CCNL 1992-2013 (https://zenodo.org/record/6644980#.Y9MQ2uzMLpA)
setwd("/Users/annadrake/Downloads/ALAN/")
##Crop and reproject and align
Lights<-list.files(pattern=".tif",all.files=TRUE, full.names=FALSE)

Light<-list()
for (j in 1:length(Light)){
  Light[[j]]<-rast(Lights[[j]])
  Light[[j]]<-terra::project(Light[[j]],Format,mask=TRUE,align=TRUE,method="near")  
  Light[[j]]<-crop(Light[[j]],e) 
  writeRaster(Light[[j]],paste(substr(Lights[j],1,4),substr(Lights[j],10,14),".tif",sep=""))
  }

#### MODIS ####
setwd("/Users/annadrake/Downloads/MODIS/")
##Crop and reproject and align
ModList<-list.files(pattern=".tif",all.files=TRUE, full.names=FALSE)
MODIS<-list()

for (j in 1:length(ModList)){
  MODIS[[j]]<-rast(ModList[[j]])
  MODIS[[j]]<-crop(MODIS[[j]],e) 
  writeRaster(MODIS[[j]],paste(ModList[j],"2.tiff",sep=""))
}
