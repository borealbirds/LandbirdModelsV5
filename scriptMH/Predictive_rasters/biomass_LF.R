################################################3
#  BAM NAM 4.1 - Producing Prediction Rasters for LANDFIRE biomass dat (height, crown closure, biomass) 
#  The script download source data from website
#  Temporal raster are saved on disk and transfered to Google SharedDrive
#################################################

library(googledrive)
library(terra)
library(dplyr)
#library(raster)
#library(data.table)

########################
#PARAM 
########################
#Extend timeout for download
options(timeout=400)

# Set extent 
rast1k <- rast(nrows=4527, ncols=7300, xmin=-4100000, xmax=3200000, ymin=1673000, ymax=6200000, crs = "EPSG:5072")

## Access urls and extract yearlist
setwd("E:/MelinaStuff/BAM/NationalModelv4.1")

webData <- read.csv("./Scripts/webDataAccess.csv", fileEncoding="UTF-8-BOM")
landfire <- subset(webData, grepl("LF",webData$Dataset))
datalist <- sort(unique(landfire$Dataset))

dwd.folder <- "./CovariateRasters/Landfire"
if (!file.exists(dwd.folder)) {
  dir.create(dwd.folder, showWarnings = FALSE)
}

# output folder
out.folder <- "./PredictionRasters/LF"
if (!file.exists(out.folder)) {
  dir.create(out.folder, showWarnings = FALSE)
}
########################
### LANDFIRE -  
########################
#DOWNLOAD
for (d in datalist){
  data <- landfire %>%
    dplyr::filter(Dataset ==d)
  yearlist <- unique(data$year)
  for(yr in yearlist){
    url <- data$url[data$year==yr]
    download.file(url, file.path(dwd.folder,paste0(d, ".zip")), method = "auto", mode = "wb")
    unzip.folder <- file.path(dwd.folder, yr)
    if (!file.exists(unzip.folder)) {
      dir.create(unzip.folder, showWarnings = FALSE)
    }
    lzip <- unzip(file.path(dwd.folder,paste0(d, ".zip")), list = TRUE)
    i <- grep(".tif$", lzip$Name )
    p<- unzip(file.path(dwd.folder,paste0(d, ".zip")), files = lzip$Name[i], exdir = unzip.folder)
    cc<-rast(p)
    writeRaster(cc, filename = file.path(unzip.folder, basename(lzip$Name[i])), overwrite=TRUE)
    #remove unecessary
    file.remove(file.path(dwd.folder, paste0(d, ".zip")))
    ls_rm <-list.dirs(unzip.folder)
    unlink(ls_rm[2:3],recursive=TRUE)
  }
}

# Resample
rast.dir <- list.dirs("./CovariateRasters/Landfire", recursive = FALSE, full.names = TRUE)
for (y in rast.dir){
  # Biomass (  0 - 40 = 0 - 0.4 kg m-3 45 = thematic class of all values > 0.4 meters AND hardwoods -9999 = NoData)
  rast_f <- list.files(y, pattern="cbd|CBD", full.names = TRUE) # Biomass
  biomass_r <- lapply(rast_f, terra::rast) #Create SpatRast
  biomass_class <- lapply(biomass_r, function(x){classify(x, cbind(-9999, NA))}) #Reclassify: value -9999 = deciduous: no value
  biomass_pj <- lapply(biomass_class, function(x) {project(x, rast1k, res = 1000, method = "bilinear", align = TRUE) }) # Reproject
  biomass_sprc <- terra::sprc(biomass_pj)  #Create SpatRastCollection
  biomass_mosaic<- terra::mosaic(biomass_sprc, fun = "max")  # Mosaic rast
  biomass <- terra::crop(biomass_mosaic, rast1k, extend = TRUE) # Extend to NatMod extent
  biomass5k <- focal(biomass_mosaic,w=matrix(1,5,5),fun=mean,na.rm=TRUE) # get 5k average at 1k resolution
  biomass5x5 <- terra::crop(biomass5k, rast1k, extend = TRUE) # Extend to NatMod extent
  #save output
  writeRaster(biomass, filename=file.path(out.folder,paste0("LFbiomass1km_", as.character(basename(y)),".tif")), overwrite=TRUE)
  writeRaster(biomass5x5, filename=file.path(out.folder,paste0("LFbiomass5x5_", as.character(basename(y)),".tif")), overwrite=TRUE)

  ##########################
  # Height
  rast_f <- list.files(y, pattern="ch|CH", full.names = TRUE)
  ht_r <- lapply(rast_f, terra::rast) #Create SpatRast
  ht_class <- lapply(ht_r, function(x){classify(x, cbind(-9999, NA))}) #Reclassify: value -9999 = deciduous: no value
  ht_pj <- lapply(ht_class, function(x) {project(x, rast1k, res = 1000, method = "bilinear", align = TRUE) }) # Reproject
  ht_sprc <- terra::sprc(ht_pj)  #Create SpatRastCollection
  ht_mosaic<- terra::mosaic(ht_sprc, fun = "max")  # Mosaic rast
  ht <- terra::crop(ht_mosaic, rast1k, extend = TRUE) # Extend to NatMod extent
  ht5k <- focal(ht_mosaic,w=matrix(1,5,5),fun=mean,na.rm=TRUE) # get 5k average at 1k resolution
  ht5x5 <- terra::crop(ht5k, rast1k, extend = TRUE) # Extend to NatMod extent

  # Heightcv
  ht_class <- lapply(ht_r, function(x){classify(x, cbind(-9999, NA))})
  ht_pj <- lapply(ht_class, function(x) {project(x, rast1k, method = "bilinear", align = TRUE) }) # Reproject
  ht_sprc <- terra::sprc(ht_pj)  #Create SpatRastCollection
  ht_mosaic<- terra::mosaic(ht_sprc, fun = "max")  # Mosaic rast
  htcv <- aggregate(ht_mosaic, fact=33, fun="sd") #coef var
  htcv_pj <- terra::project(htcv, rast1k)
  htcv_cp <- terra::crop(htcv_pj, rast1k, extend = TRUE) # Extend to NatMod extent
  htcv5k<-terra::focal(htcv_pj,w=matrix(1,5,5),fun=mean,na.rm=TRUE) #get 5k average at 1k resolution
  htcv5k_cp <- terra::crop(htcv5k, rast1k, extend = TRUE) # Extend to NatMod extent
  
  writeRaster(htcv_cp, filename=file.path(out.folder,paste0("LFheigthcv1km_", as.character(basename(y)),".tif")), overwrite=TRUE)
  writeRaster(htcv5k_cp, filename=file.path(out.folder,paste0("LFheigthcv5k_", as.character(basename(y)),".tif")), overwrite=TRUE)
  
  #save output
  writeRaster(ht, filename=file.path(out.folder,paste0("LFheigth1km_", as.character(basename(y)),".tif")), overwrite=TRUE)
  writeRaster(ht5x5, filename=file.path(out.folder,paste0("LFheigth5x5_", as.character(basename(y)),".tif")), overwrite=TRUE)
  writeRaster(htcv, filename=file.path(out.folder,paste0("LFheigthcv1km_", as.character(basename(y)),".tif")), overwrite=TRUE)
  writeRaster(htcv5k, filename=file.path(out.folder,paste0("LFheigthcv5k_", as.character(basename(y)),".tif")), overwrite=TRUE)
  
   ##########################
  # Crown Closure
  rast_f <- list.files(y, pattern="cc|CC", full.names = TRUE)
  cc_r <- lapply(rast_f, terra::rast) #Create SpatRast
  cc_class <- lapply(cc_r, function(x){classify(x, cbind(-9999, NA))}) #Reclassify: value -9999 = deciduous: no value
  cc_pj <- lapply(cc_class, function(x) {project(x, rast1k, res = 1000, method = "bilinear", align = TRUE) }) # Reproject
  cc_sprc <- terra::sprc(cc_pj)  #Create SpatRastCollection
  cc_mosaic<- terra::mosaic(cc_sprc, fun = "max")  # Mosaic rast
  cc <- terra::crop(cc_mosaic, rast1k, extend = TRUE) # Extend to NatMod extent
  cc5k <- focal(cc_mosaic,w=matrix(1,5,5),fun=mean,na.rm=TRUE) # get 5k average at 1k resolution
  cc5x5 <- terra::crop(cc5k, rast1k, extend = TRUE) # Extend to NatMod extent
  #save output
  writeRaster(cc, filename=file.path(out.folder,paste0("LFcrownclosure1km_", as.character(basename(y)),".tif")), overwrite=TRUE)
  writeRaster(cc5x5, filename=file.path(out.folder,paste0("LFcrownclosure5x5_", as.character(basename(y)),".tif")), overwrite=TRUE)
}
  
  
  
