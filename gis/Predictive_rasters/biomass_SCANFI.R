################################################3
#  BAM NAM 4.1 - Producing Prediction Rasters for SCANFI biomass data 
#  The script download source data from the sharedDrive
#  Temporal raster are saved on disk and transferred to Google SharedDrive
#################################################

library(googledrive)
library(terra)


########################
### PARAMS
########################
# workdir
setwd("E:/MelinaStuff/BAM/NationalModelv5.0")

# Set extent 
rast1k <- rast(nrows=4527, ncols=7300, xmin=-4100000, xmax=3200000, ymin=1673000, ymax=6200000, crs = "EPSG:5072")

# output folder
out_f <- "./PredictionRasters/SCANFI2"
if (!file.exists(out_f)) {
  dir.create(out_f, showWarnings = FALSE)
}

## Access urls and extract yrlist
webData <- read.csv("./Scripts/webDataAccess.csv", fileEncoding="UTF-8-BOM")
SCANFI_lnd <- subset(webData, Dataset == "Biomass_SCANFI")
SCANFI_parentid <- SCANFI_lnd$url

########################
### SCANFI SOURCE
########################
#Extract GoogleDrive id to store output
SCANFI_list <- c("SCANFI_1985", "SCANFI_1990", "SCANFI_1995", "SCANFI_2000",
                 "SCANFI_2005", "SCANFI_2010", "SCANFI_2015", "SCANFI_2020")

dr <- drive_ls(as_id(SCANFI_parentid), pattern = "SCANFI")

##########################
##BIOMASS
#Download biomass from GDrive
for (j in 1:length(SCANFI_list)){
  dr_subf <- drive_ls(as_id(SCANFI_parentid), pattern = SCANFI_list[j])
  dr_bio <- drive_ls(as_id(dr_subf), pattern = "Biomass")
  for (i in 1:nrow(dr_bio)){
    if(endsWith(dr_bio$name[i], "tif"))
      drive_download(as_id(dr_bio$id[i]), path = file.path("./CovariateRasters/SCANFI", dr_bio$name[i]))
  }
}

#Generate predictive raster
biomass_list <- list.files(path = "./CovariateRasters/SCANFI/", pattern = "Biomass*", full.names = TRUE)
for (j in 1:length(biomass_list)){
  outfile <- basename(biomass_list[j])
  biomass<-rast(biomass_list[j])
  biomass_pj <- terra::project(biomass, rast1k)
  biomass5k<-terra::focal(biomass_pj,w=matrix(1,5,5),fun=mean,na.rm=TRUE) #get 5k average at 1k resolution
  
  biomass_cp<-crop(biomass_pj, ext(rast1k))
  biomass5k_cp<-crop(biomass5k, ext(rast1k))
  
  writeRaster(biomass_cp, file.path(out_f, paste0("SCANFIbiomass1km_", substr(outfile, 22, 25),".tif")), overwrite=TRUE)
  writeRaster(biomass5k_cp, file.path(out_f, paste0("SCANFIbiomass5x5_", substr(outfile, 22, 25),".tif")), overwrite=TRUE)
}

##########################
##HEIGHT
# Google drive batch download
for (j in 1:length(SCANFI_list)){
  dr_subf <- drive_ls(as_id(SCANFI_parentid), pattern = SCANFI_list[j])
  dr_bio <- drive_ls(as_id(dr_subf), pattern = "height")
  for (i in 1:nrow(dr_bio)){
    if(endsWith(dr_bio$name[i], "tif"))
      drive_download(as_id(dr_bio$id[i]), path = file.path("./CovariateRasters/SCANFI", dr_bio$name[i]))
  }
}

# batch process
#cv
cv <- function(x){
  y <- na.omit(sample(x, size=10, replace =F))
  sd(y) / mean(y)
}

height_list <- list.files(path = "./CovariateRasters/SCANFI/", pattern = "height*", full.names = TRUE)
for (j in 1:length(height_list)){
  outfile <- basename(height_list[j])
  height<-rast(height_list[j])
  height_pj <- terra::project(height, rast1k)
  height5k<-terra::focal(height_pj,w=matrix(1,5,5),fun=mean,na.rm=TRUE) #get 5k average at 1k resolution
  heightcv <- aggregate(height, fact=33, fun=cv)
  heightcv_pj <- terra::project(heightcv, rast1k, res = 1000, method = "bilinear", align = TRUE)
  heightcv5k<-terra::focal(heightcv_pj,w=matrix(1,5,5),fun=mean,na.rm=TRUE) #get 5k average at 1k resolution
  
  height_cp<-crop(height_pj, ext(rast1k))
  height5k_cp<-crop(height5k, ext(rast1k))
  heightcv_cp<-crop(heightcv_pj, ext(rast1k))
  heightcv5k_cp<-crop(heightcv5k, ext(rast1k))
  
  writeRaster(height_cp, file.path(out_f, paste0("SCANFIheight1km_", substr(outfile, 21, 24),".tif")), overwrite=TRUE)
  writeRaster(height5k_cp, file.path(out_f, paste0("SCANFIheight5x5_", substr(outfile, 21, 24),".tif")), overwrite=TRUE)
  writeRaster(heightcv_cp, file.path(out_f, paste0("SCANFIheightcv1km_", substr(outfile, 21, 24),".tif")), overwrite=TRUE)
  writeRaster(heightcv5k_cp, file.path(out_f, paste0("SCANFIheightcv5x5_", substr(outfile, 21, 24),".tif")), overwrite=TRUE)
}


##########################
##CROWN CLOSURE
# Google drive batch download
for (j in 1:length(SCANFI_list)){
  dr_subf <- drive_ls(as_id(SCANFI_parentid), pattern = SCANFI_list[j])
  dr_bio <- drive_ls(as_id(dr_subf), pattern = "closure")
  for (i in 1:nrow(dr_bio)){
    if(endsWith(dr_bio$name[i], "tif"))
      drive_download(as_id(dr_bio$id[i]), path = file.path("./CovariateRasters/SCANFI", dr_bio$name[i]))
  }
}

# batch process
cc_list <- list.files(path = "./CovariateRasters/SCANFI/", pattern = "closure*", full.names = TRUE)
for (j in 1:length(cc_list)){
  outfile <- basename(cc_list[j])
  cc<-rast(cc_list[j])
  cc_pj <- terra::project(cc, rast1k)
  cc5k<-terra::focal(cc_pj,w=matrix(1,5,5),fun=mean,na.rm=TRUE) #get 5k average at 1k resolution
  
  cc_cp<-crop(cc_pj, e)
  cc5k_cp<-crop(cc5k, e)
  
  writeRaster(cc_cp, file.path(out_f, paste0("SCANFIclosure1km_", substr(outfile, 22, 25),".tif")), overwrite=TRUE)
  writeRaster(cc5k_cp, file.path(out_f, paste0("SCANFIclosure5x5_", substr(outfile, 22, 25),".tif")), overwrite=TRUE)
}


##########################
##   TREE SPECIES
# Google drive batch download
for (j in 1:length(SCANFI_list)){
  dr_subf <- drive_ls(as_id(SCANFI_parentid), pattern = SCANFI_list[j])
  dr_bio <- drive_ls(as_id(dr_subf), pattern = "sps")
  for (i in 1:nrow(dr_bio)){
    if(endsWith(dr_bio$name[i], "tif"))
      drive_download(as_id(dr_bio$id[i]), path = file.path("./CovariateRasters/SCANFI", dr_bio$name[i]))
  }
}

# batch process
cc_list <- list.files(path = "./CovariateRasters/SCANFI/", pattern = "sps*", full.names = TRUE)
for (j in 1:length(cc_list)){
  outfile <- basename(cc_list[j])
  cc<-rast(cc_list[j])
  cc_pj <- terra::project(cc, rast1k)
  cc5k<-terra::focal(cc_pj,w=matrix(1,5,5),fun=mean,na.rm=TRUE) #get 5k average at 1k resolution
  
  cc_cp<-crop(cc_pj, e)
  cc5k_cp<-crop(cc5k, e)
  
  if(strsplit(outfile ,"_")[[1]][5] %in% c("1985","1990","1995","2000","2005","2010","2015","2020")){
    writeRaster(cc_cp, file.path(out_f, paste0("SCANFI", strsplit(outfile ,"_")[[1]][3], "1km", strsplit(outfile ,"_")[[1]][5],".tif")), overwrite=TRUE)
    writeRaster(cc5k_cp, file.path(out_f, paste0("SCANFI", strsplit(outfile ,"_")[[1]][3], "5x5", strsplit(outfile ,"_")[[1]][5],".tif")), overwrite=TRUE)
  }else{
    writeRaster(cc_cp, file.path(out_f, paste0("SCANFI", strsplit(outfile ,"_")[[1]][3], "1km", strsplit(outfile ,"_")[[1]][6],".tif")), overwrite=TRUE)
    writeRaster(cc5k_cp, file.path(out_f, paste0("SCANFI", strsplit(outfile ,"_")[[1]][3], "5x5", strsplit(outfile ,"_")[[1]][6],".tif")), overwrite=TRUE)
  }
}
