################################################3
#  BAM NM 5.0 - Producing Prediction Rasters for Sentinel canopy height  data for the year 2020
#  The script download from GEE, reproject and rescale to 1km
#  Ouptput raster is saved on disk and transfered to Google SharedDrive
#################################################

#install.packages('terra', repos='https://rspatial.r-universe.dev')
library(terra)
library(sp)
library(rgee)
library(reticulate)
library(googledrive)

# Initialize Earth Engine
ee_check()
ee_Authenticate()
ee_Initialize(user = "houle.melina@gmail.com")
googledrive::drive_auth("houle.melina@gmail.com")



# workdir
setwd("E:/MelinaStuff/BAM/NationalModelv5.0")

# Set extent 
rast1k <- rast(nrows=4527, ncols=7300, xmin=-4100000, xmax=3200000, ymin=1673000, ymax=6200000, crs = "EPSG:5072")

# Set output and download folder
out_f <- "./PredictionRasters/ETH2"
out_f <- "./temp/ETH2"

if (!file.exists(out_f)) {
  dir.create(out_f, showWarnings = FALSE)
}


## Access urls 
webData <- read.csv("./Scripts/webDataAccess.csv", fileEncoding="UTF-8-BOM")
ETH_lnd <- subset(webData, Dataset == "Sentinel canopy height")
img_name <- ETH_lnd$url

##########################################################
##Set function
##########################################################
##Set function
ETHheight <- function(x) {
  ccc <- ee$Image(x)$ # read the data
    select('b1')#
  return(ccc)
}

##########################################################
# Extract image on GEE per year
setwd(out_f)
imageextract <- ETHheight(img_name)
geomPoly <- ee$Geometry$BBox(-180, 38, -45, 70);

temp <- ee_as_rast(
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
r_ETHheight <- project(raster_mosaic, rast1k, res = 1000, method = "bilinear", align = TRUE) # Reproject
ETHheight <- crop(r_ETHheight, ext(rast1k), extend = TRUE) # Crop to NatMod extent
#5x5
ETHfocal5k<-terra::focal(r_ETHheight,w=matrix(1,5,5),fun=mean,na.rm=TRUE) #get 5k average at 1k resolution
ETH5k <- terra::crop(ETHfocal5k, rast1k, extend = TRUE) # Extend to NatMod extent

#cv
cv <- function(x){
  y <- na.omit(sample(x, size=10, replace =F))
  sd(y) / mean(y)
}

r_ETHheightcv <- aggregate(r_ETHheight, fact=10, fun=cv) #coefficient of variation
r_ETHheightcv1k <- project(r_ETHheightcv, rast1k, res = 1000, method = "bilinear", align = TRUE) # Reproject
ETHheightcv1k <- terra::crop(r_ETHheightcv1k, rast1k, extend = TRUE) # Extend to NatMod extent
r_ETHheightcv5k<-terra::focal(r_ETHheightcv1k,w=matrix(1,5,5),fun=mean,na.rm=TRUE) #get 5k average at 1k resolution
ETHheightcv5k <- terra::crop(r_ETHheightcv5k, rast1k, extend = TRUE) # Extend to NatMod extent


#write output
writeRaster(ETHheight, filename="./ETHheight.tif", filetype="GTiff", overwrite=TRUE)
writeRaster(ETH5k, filename="./ETHheight5k.tif", filetype="GTiff", overwrite=TRUE)
writeRaster(ETHheightcv1k, filename="./ETHheightcv.tif", filetype="GTiff", overwrite=TRUE)
writeRaster(ETHheightcv5k, filename="./ETHheightcv5k.tif", filetype="GTiff", overwrite=TRUE)
