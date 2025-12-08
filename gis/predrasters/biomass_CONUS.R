################################################3
#  BAM NM 5.0 - Producing Prediction Rasters for CONUS biomass tree cover 
#  The script download from GEE, reproject and rescale to 100m and 1km
#  Temporal raster are saved on disk and transfered to Google SharedDrive
#################################################
#####################################################
##   CONUS CROWN CLOSURE
##
#####################################################
library(googledrive)
library(terra)
library(reticulate)
library(rgee)

##########################################################
#PARAM 
Sys.setenv(HOME="C:/Users/mehou10/Documents")
setwd("E:/MelinaStuff/BAM/NationalModelv5.0")


# Set extent 
rast1k <- rast(nrows=4527, ncols=7300, xmin=-4100000, xmax=3200000, ymin=1673000, ymax=6200000, crs = "EPSG:5072")

#use_python("C:/Users/mehou10/AppData/Local/r-miniconda/envs/r-reticulate")
ee_Initialize(user = "houle.melina@gmail.com", drive = TRUE)
ee_check()

## Access urls and extract yrlist
webData <- read.csv("./Scripts/webDataAccess.csv", fileEncoding="UTF-8-BOM")
CONUS <- subset(webData, Dataset == "Biomass_CONUS")
yrlist <- seq(as.numeric(strsplit(CONUS$year, "[:]")[[1]][1]), as.numeric(strsplit(CONUS$year, "[:]")[[1]][2]),by=1)
img_name <- CONUS$url

out_f <- "./PredictionRasters/CONUS"
if (!file.exists(out_f)) {
  dir.create(out_f, showWarnings = TRUE, recursive = TRUE)
}
##########################################################
##Set function
conus_cc <- function(yr) {
  ccc <- ee$ImageCollection(img_name)$ # read the data
    select('TRE')$ # select landcover
    filter(ee$Filter$calendarRange(yr, yr, "year"))$
    toBands();# filter the year
  return(ccc)
}

##########################################################

# Extract image on GEE per year
for (i in 1:length(yrlist)) {
  imageextract <- conus_cc(as.integer(yrlist[i]))
  temp <- ee_as_raster(
    image = imageextract,
    container = "RAP_TCC",
    dsn = paste0("RAP_TCC1km_",as.character(yrlist[i])),
    scale = 1000,
    crs = 'EPSG:5072',
    via = "drive",
    maxPixels= 4e10
  )
  # Convert to SPatRast (terra format)
  spat_rast <- rast(temp)
  # Reproject
  spat_rast.pj <- project(spat_rast, rast1km, res = 1000, method = "near", align = TRUE)
  # Crop to NatMod extent
  rast.cp <- crop(spat_rast.pj, ext(rast1km), extend = TRUE)
  
  #write output
  writeRaster(rast.cp, filename=file.path(out_f,paste0("RAP_TCC1km_",as.character(yrlist[i]), ".tif")), filetype="GTiff", overwrite=TRUE)
  file.remove(rastlist)
}

#get 5k average at 1k resolution
rastlist <- list.files(path = out_f, pattern="RAP_TCC1km_", 
                       all.files=TRUE, full.names=FALSE)
rasters_1km <-lapply(file.path(out_f, rastlist), rast)
rasters5k<-lapply(rasters_1km, function(x) {focal(x,w=matrix(1,5,5),fun=mean,na.rm=TRUE)})
for (i in 1:length(rastlist)) {names(rasters5k[[i]]) <- paste0("RAP_TCC5x5_", substr(rastlist[i],12,15), ".tif")}
lapply(rasters5k, function(i) {
  writeRaster(i, filename=file.path(out_f,names(i)), filetype="GTiff", overwrite = TRUE)
})







###################################################################
##### Resolution  100m
# Extract image on GEE per year
for (i in 1:length(yrlist)) {
  imageextract <- conus_cc(as.integer(yrlist[i]))
  temp <- ee_as_raster(
    image = imageextract,
    container = "RAP_TCC",
    dsn = paste0("RAP_TCC1km_",as.character(yrlist[i])),
    scale = 1000,
    crs = 'EPSG:5072',
    via = "drive",
    maxPixels= 4e10
  )
  #List downloaded images fo a specific year
  rastlist <- list.files(path = "./", pattern=paste0("RAP_TCC100m_", yrlist[i]), 
                         all.files=TRUE, full.names=FALSE)
  # Load images as SpatRaster
  flf <-lapply(rastlist, rast)
  # Create SpatRastCollection
  raster_sprc <- terra::sprc(flf)
  # Mosaic all SpatRast together
  raster_mosaic<- terra::mosaic(raster_sprc, fun = "mean") 
  # Reproject
  r_b1pj <- project(raster_mosaic, rast100, res = 100, method = "near", align = TRUE)
  # Crop to NatMod extent
  rast <- crop(r_b1pj, ext(rast100), extend = TRUE)
  
  #write output
  writeRaster(rast, filename=paste0("RAP_TCC100m_",as.character(yrlist[i]), ".tif"), filetype="GTiff", overwrite=TRUE)
  file.remove(rastlist)
}