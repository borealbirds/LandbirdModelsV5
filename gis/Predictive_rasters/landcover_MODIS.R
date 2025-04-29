################################################3
#  BAM NM 5.0 - Producing Prediction Rasters for MODIS LCC 2001-2020
#  The script download from GEE, reproject and rescale to 1km with a 5x5 moving window
#  Ouptput raster is saved on disk and transfered to Google SharedDrive
#################################################

#install.packages('terra', repos='https://rspatial.r-universe.dev')
library(terra)
library(sp)
library(rgee)
library(reticulate)
library(googledrive)
library(future)

# Initialize Earth Engine
ee_check()
ee_Initialize(user = "houle.melina@gmail.com")

drive_auth(email= "houle.melina@gmail.com")

## workdir
setwd("E:/MelinaStuff/BAM/NationalModelv5.0")

# Set extent 
rast1k <- rast(nrows=4527, ncols=7300, xmin=-4100000, xmax=3200000, ymin=1673000, ymax=6200000, crs = "EPSG:5072")

rast1k_wgs84 <- project(rast1k, "EPSG:4326")
ext <- ext(rast1k_wgs84)

# Define GEE geometry using the projected extent
region_geom <- ee$Geometry$Rectangle(
  coords = c(ext[1], ext[3], ext[2], ext[4]),  # xmin, ymin, xmax, ymax
  proj = "EPSG:4326",
  geodesic = FALSE
)

# Set output and download folder
out_f <- "./PredictionRasters/MODIS"
if (!file.exists(out_f)) {
  dir.create(out_f, showWarnings = FALSE)
}

## Access urls 
webData <- read.csv("./Scripts/webDataAccess.csv", fileEncoding="UTF-8-BOM")
MODIS_lnd <- subset(webData, Dataset == "MODISlcc")
img_name <- MODIS_lnd$url
#yrlist <- seq(as.numeric(strsplit(MODIS_lnd$year, "[:]")[[1]][1]), as.numeric(strsplit(MODIS_lnd$year, "[:]")[[1]][2]),by=1)

yrlist <- eval(parse(text = MODIS_lnd$year[1]))

##########################################################
##Set function
##########################################################
##Set function
MODIS500 <- function(yr) {
  modis_img <- ee$ImageCollection("MODIS/006/MCD12Q1")$
    filterDate(paste0(yr, "-01-01"), paste0(yr, "-12-31"))$
    first()$
    select("LC_Type1")$
    reproject(crs = "EPSG:4326", scale = 500)
  
  img_rast <- ee_as_rast(
    image = modis_img,
    region = region_geom,
    scale = 500,
    via = "drive",
    dsn = paste0("MODIS500m_epsg4326_", yr, ".tif"),
    container = "rgee_backup",
    maxPixels = 1e9,
    lazy = FALSE,
    public = TRUE,
    timePrefix = FALSE,
    quiet = FALSE
  )
  
  return(img_rast)
}
# Run function to extract
for (i in 1:length(yrlist)) {
  MODISextract <- MODIS500(yrlist[i])
  MODISlist <- list.files(path = "./", pattern=paste0("MODIS500m_epsg4326_", as.character(yrlist[i])), 
                         all.files=TRUE, full.names=TRUE)
  flf <-lapply(MODISlist, rast) # Load images as SpatRaster
  raster_sprc <- terra::sprc(flf) # Create SpatRastCollection
  raster_mosaic<- terra::mosaic(raster_sprc, fun = "max")  # Mosaic all SpatRast together
  MODIS_pj <- project(raster_mosaic, rast1k, res = 1000, method = "near", align = TRUE) # Reproject
  MODIS_crop <- crop(MODIS_pj, rast1k)
  writeRaster(MODIS_crop, filename=file.path(out_f, paste0("MODISLCC_1km_",as.character(yrlist[i]), ".tif")), filetype="GTiff", overwrite=TRUE)
  MODISfocal5k<-terra::focal(MODIS_pj, w=matrix(1,5,5),fun=modal, na.rm=TRUE) #get 5k average at 1k resolution
  MODISlcc5k <-crop(MODISfocal5k, rast1k, extend = TRUE)
  writeRaster(MODISlcc5k, filename=file.path(out_f, paste0("MODISLCC_5x5_",as.character(yrlist[i]), ".tif")), filetype="GTiff", overwrite=TRUE)
}



#------------
#List downloaded images fo a specific year


flf <-lapply(rastlist, rast) # Load images as SpatRaster
raster_sprc <- terra::sprc(flf) # Create SpatRastCollection
raster_mosaic<- terra::mosaic(raster_sprc, fun = "mean")  # Mosaic all SpatRast together
r_abgpj <- project(raster_mosaic, rast1km, res = 1000, method = "bilinear", align = TRUE) # Reproject
rast <- crop(r_abgpj, ext(rast1km), extend = TRUE) # Crop to NatMod extent

ONRLfocal5k<-terra::focal(r_abgpj,w=matrix(1,5,5),fun=mean,na.rm=TRUE) #get 5k average at 1k resolution
ONRL5k <- terra::crop(ONRLfocal5k, rast1k, extend = TRUE) # Extend to NatMod extent

#write output
writeRaster(rast, filename=file.path(out_f,"ONRLbiomass.tif"), filetype="GTiff", overwrite=TRUE)
writeRaster(ONRL5k, filename=file.path(out_f,"ONRLbiomass5k.tif"), filetype="GTiff", overwrite=TRUE)










##########################################################
# Extract image on GEE per year
setwd(out_f)
imageextract <- ETHheight(img_name)
geomPoly <- ee$Geometry$BBox(-180, 38, -45, 70);

temp <- ee_as_raster(
    image = imageextract,
    region = geomPoly,
    container = "ETH",
    dsn = "ETH",
    scale = 100,
    #crs = 'EPSG:5072',
    via = "drive",
    maxPixels= 6e12
  )

#List downloaded images fo a specific year
rastlist <- list.files(path = "./", pattern="ETH", 
                       all.files=TRUE, full.names=FALSE)

flf <-lapply(rastlist, rast) # Load images as SpatRaster
raster_sprc <- terra::sprc(flf) # Create SpatRastCollection
raster_mosaic<- terra::mosaic(raster_sprc, fun = "mean")  # Mosaic all SpatRast together
r_ETHheight <- project(raster_mosaic, rast1km, res = 1000, method = "bilinear", align = TRUE) # Reproject
ETHheight <- crop(r_ETHheight, ext(rast1km), extend = TRUE) # Crop to NatMod extent
#5x5
ETHfocal5k<-terra::focal(r_ETHheight,w=matrix(1,5,5),fun=mean,na.rm=TRUE) #get 5k average at 1k resolution
ETH5k <- terra::crop(ETHfocal5k, rast1k, extend = TRUE) # Extend to NatMod extent

#cv
r_ETHheight <- project(raster_mosaic, rast1km, method = "bilinear", align = TRUE) # Reproject
r_ETHheightcv <- aggregate(r_ETHheight, fact=10, fun="sd") #coef var
r_ETHheightcv1k <- project(r_ETHheightcv, rast1km, res = 1000, method = "bilinear", align = TRUE) # Reproject
ETHheightcv1k <- terra::crop(r_ETHheightcv1k, rast1k, extend = TRUE) # Extend to NatMod extent
r_ETHheightcv5k<-terra::focal(r_ETHheightcv1k,w=matrix(1,5,5),fun=mean,na.rm=TRUE) #get 5k average at 1k resolution
ETHheightcv5k <- terra::crop(r_ETHheightcv5k, rast1k, extend = TRUE) # Extend to NatMod extent

#write output
writeRaster(ETHheight, filename="./ETHheight.tif", filetype="GTiff", overwrite=TRUE)
writeRaster(ETH5k, filename="./ETHheight5k.tif", filetype="GTiff", overwrite=TRUE)
writeRaster(ETHheightcv1k, filename="./ETHheightcv.tif", filetype="GTiff", overwrite=TRUE)
writeRaster(ETHheightcv5k, filename="./ETHheightcv5k.tif", filetype="GTiff", overwrite=TRUE)

