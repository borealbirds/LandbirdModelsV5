######################################################
# Title: Recalculate CV for Landfire and ETH/Sentinel using equivalent scale to covariate extraction
# Author: Anna Drake
# Date: June 2024
########################################################
# Notes - this resolves persistant mismatch/bias between prediction values extracted at 1km 
# and covariate values. The scale of the calculation between the two outputs is brought into alignment
# and then the prediction layers are set to the correct resolution.

require(raster)
require(terra)
require(dplyr)
require(sf)

#Set root path for data on google drive ----
root<-"G:/Shared drives/BAM_NationalModels5" #PC link

# Template layer ----
r<-terra::rast(file.path(root,"PredictionRasters","Topography","mTPI_1km.tif"))

#1. LF height_cv ----
for (i in c(2001,2012,2014,2022)){

#2. ALASKA ----
# Import 30m res ALASKA  --------------
path<-file.path(root,"CovariateRasters","Biomass","Landfire","sourceData","AK",i)
file<-list.files(path,pattern="ch|CH", full.names = TRUE)
height<-terra::rast(file)
ht<-classify(height, cbind(-9999, NA)) #Reclassify: value -9999 = deciduous: no value

# ~200m radius = 350m2 ----
mean200<-aggregate(ht, fact=11, fun="mean", na.rm=T)
mean200[mean200<0.1]<-0
sd200<-aggregate(ht, fact=11, fun="sd", na.rm=T)
LFcv200<-sd200/mean200 
LFcv200[!is.finite(LFcv200)]<-NA
LFcv1km_AK<-project(LFcv200,r, method="bilinear", align =T)

# ~2000m radius = 3500m2 ----
mean2000<-aggregate(ht, fact=117, fun="mean", na.rm=T)
mean2000[mean2000<0.1]<-0
sd2000<-aggregate(ht, fact=117, fun="sd", na.rm=T)
LFcv2000<-sd2000/mean2000 
LFcv2000[!is.finite(LFcv2000)]<-NA
LFcv5km_AK<-disagg(LFcv2000,fact=5, method="bilinear")%>% project(.,r, method="bilinear", align =T)

#3. CONUS ----
# Import 30m res file CONUS ----   
path<-file.path(root,"CovariateRasters","Biomass","Landfire","sourceData","CONUS",i)
file<-list.files(path,pattern="ch|CH", full.names = TRUE)
height<-terra::rast(file)
ht<-classify(height, cbind(-9999, NA)) #Reclassify: value -9999 = deciduous: no value

# ~200m radius = 350m2 ----
mean200<-aggregate(ht, fact=11, fun="mean", na.rm=T)
mean200[mean200<0.1]<-0
sd200<-aggregate(ht, fact=11, fun="sd", na.rm=T)
LFcv200<-sd200/mean200  
LFcv200[!is.finite(LFcv200)]<-NA
LFcv1km_CONUS<-project(LFcv200,r, method="bilinear", align =T)

# ~2000m radius = 3500m2 ----
mean2000<-aggregate(ht, fact=117, fun="mean", na.rm=T)
mean2000[mean2000<0.1]<-0
sd2000<-aggregate(ht, fact=117, fun="sd", na.rm=T)
LFcv2000<-sd2000/mean2000 
LFcv2000[!is.finite(LFcv2000)]<-NA
LFcv5km_CONUS<-disagg(LFcv2000,fact=5, method="bilinear")%>% project(.,r, method="bilinear", align =T)

### Save final tiff files ----
Final1k<-merge(LFcv1km_AK,LFcv1km_CONUS) %>% crop(.,r) 
writeRaster(Final1k,paste("LF_cv_recal",i,".tif",sep=""))

Final5x5<-merge(LFcv5km_CONUS,LFcv5km_AK) %>% crop(.,r) 
writeRaster(Final5x5,paste("LF_cv_5x5_recal",i,".tif",sep=""))

} # all years

#4. Sentinel 2020 10m exported as the 30m mean -----
# Exported from GEE as (e.g. Prairie region):
#var region = ee.Geometry.BBox(-119,48,-92,56);
#var canopy_height = ee.Image("users/nlang/ETH_GlobalCanopyHeight_2020_10m_v1");
#Export.image.toDrive({
#  image: canopy_height,
#  description: 'Sentinal_ht',
#  scale: 30,
#  region: region,
#  maxPixels: 1e13,
#  crs: 'EPSG:5072'
#});
   
Import files ---------------
ETH_AK<-terra::rast("C:/Users/andrake/Documents/ETH/Sentinal_htAK.tif")
ETH_Prairie<-terra::rast("C:/Users/andrake/Documents/ETH/Sentinal_ht-0000000000-0000065536.tif")
ETH_Prairie2<-terra::rast("C:/Users/andrake/Documents/ETH/Sentinal_ht-0000000000-0000000000.tif")
ETH_Prairie<-merge(ETH_Prairie,ETH_Prairie2)

#5. ALASKA layer -----

#~200m radius = 350m2 -----
mean200<-aggregate(ETH_AK, fact=11, fun="mean", na.rm=T)
mean200[mean200<0.1]<-0
sd200<-aggregate(ETH_AK, fact=11, fun="sd", na.rm=T)
ETHcv200<-sd200/mean200 
ETHcv200[!is.finite(ETHcv200)]<-NA
ETHcv1km_AK<-project(ETHcv200,r, method="bilinear", align =T)

# ~2000m radius = 3500m2 -----
mean2000<-aggregate(ETH_AK, fact=117, fun="mean", na.rm=T)
mean2000[mean2000<0.1]<-0
sd2000<-aggregate(ETH_AK, fact=117, fun="sd", na.rm=T)
ETHcv2000<-sd2000/mean2000 
ETHcv2000[!is.finite(ETHcv2000)]<-NA
ETHcv5km_AK<-disagg(ETHcv2000,fact=5, method="bilinear")%>% project(.,r, method="bilinear", align =T)

#6. Prairie layer -----

#~200m radius = 350m2 ------
mean200<-aggregate(ETH_Prairie, fact=11, fun="mean", na.rm=T)
mean200[mean200<0.1]<-0
sd200<-aggregate(ETH_Prairie, fact=11, fun="sd", na.rm=T)
ETHcv200<-sd200/mean200 
ETHcv200[!is.finite(ETHcv200)]<-NA
ETHcv1km_P<-project(ETHcv200,r, method="bilinear", align =T)

#~2000m radius = 3500m2 ------
mean2000<-aggregate(ETH_Prairie, fact=117, fun="mean", na.rm=T)
mean2000[mean2000<0.1]<-0
sd2000<-aggregate(ETH_Prairie, fact=117, fun="sd", na.rm=T)
ETHcv2000<-sd2000/mean2000 
ETHcv2000[!is.finite(ETHcv2000)]<-NA
ETHcv5km_P<-disagg(ETHcv2000,fact=5, method="bilinear")%>% project(.,r, method="bilinear", align =T)

## Write out files -------------
ETH_Final1k<-merge(ETHcv1km_AK,ETHcv1km_P) %>% crop(.,r) %>% extend(.,r, fill=NA)
writeRaster(ETH_Final1k,"ETH_cv_1km_recal.tif", overwrite=T)

ETH_Final5k<-merge(ETHcv5km_AK,ETHcv5km_P) %>% crop(.,r) %>% extend(.,r, fill=NA)
writeRaster(ETH_Final5k,"ETH_cv_5x5_recal.tif", overwrite=T)

####################### END OF SCRIPT ####################################

