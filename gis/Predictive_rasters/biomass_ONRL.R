################################################3
#  BAM NM 5.0 - Producing Prediction Rasters for biomass ONRL carbon density data for the year 2010
#  The script download from GEE, reproject and rescale to 1km
#  Ouptput raster is saved on disk and transfered to Google SharedDrive
#################################################

#install.packages('terra', repos='https://rspatial.r-universe.dev')
library(terra)
library(sp)
library(rgee)
library(reticulate)

# Initialize Earth Engine
ee_check()
ee_Initialize()

# workdir
setwd("E:/MelinaStuff/BAM/NationalModelv5.0")

# Set extent 
rast1k <- rast(nrows=4527, ncols=7300, xmin=-4100000, xmax=3200000, ymin=1673000, ymax=6200000, crs = "EPSG:5072")

# Set output and download folder
out_f <- "./PredictionRasters/ORNL"
if (!file.exists(out_f)) {
  dir.create(out_f, showWarnings = FALSE)
}

## Access urls and extract yrlist
webData <- read.csv("./Scripts/webDataAccess.csv", fileEncoding="UTF-8-BOM")
carbon_lnd <- subset(webData, Dataset == "Biomass carbon density")
img_name <- carbon_lnd$url

##########################################################
##Set function
##########################################################
##Set function
carbonDen <- function(x) {
  ccc <- ee$ImageCollection(x)$ # read the data
    select('agb')$ # select landcover
    toBands();# filter the year
  return(ccc)
}

##########################################################
# Extract image on GEE per year
imageextract <- carbonDen(img_name)
temp <- ee_as_raster(
    image = imageextract,
    container = "ONRL",
    dsn = "ONRL_1km_",
    scale = 300,
    #crs = 'EPSG:5072',
    via = "drive",
    maxPixels= 4e10
  )
#List downloaded images fo a specific year
rastlist <- list.files(path = "./", pattern="ONRL_1km_", 
                       all.files=TRUE, full.names=FALSE)

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


