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

# Initialize Earth Engine
ee_check()
ee_Initialize()

# workdir
setwd("E:/MelinaStuff/BAM/NationalModelv4.1")

# Set extent 
rast1km <- rast(nrows=4527, ncols=7300, xmin=-4100000, xmax=3200000, ymin=1673000, ymax=6200000, crs = "EPSG:5072")

# Set output and download folder
out_f <- "./PredictionRasters/ETH"
if (!file.exists(out_f)) {
  dir.create(out_f, showWarnings = FALSE)
}
setwd(out_f)

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

#cv 5x5

#write output
writeRaster(ETHheight, filename="./ETHheight.tif", filetype="GTiff", overwrite=TRUE)
writeRaster(ETH5k, filename="./ETHheight5k.tif", filetype="GTiff", overwrite=TRUE)
writeRaster(ETHheightcv1k, filename="./ETHheightcv.tif", filetype="GTiff", overwrite=TRUE)
writeRaster(ETHheightcv5k, filename="./ETHheightcv5k.tif", filetype="GTiff", overwrite=TRUE)






# Convert to SPatRast (terra format)
spat_rast <- rast(temp)
# Reproject
spat_rast.pj <- project(spat_rast, rast1km, res = 1000, method = "bilinear", align = TRUE)
# Crop to NatMod extent
rast.cp <- crop(spat_rast.pj, ext(rast1km), extend = TRUE)
  
#write output
writeRaster(rast.cp, filename=file.path(out_f,paste0("RAP_TCC1km_",as.character(yrlist[i]), ".tif")), filetype="GTiff", overwrite=TRUE)
file.remove(rastlist)


#get 5k average at 1k resolution
rastlist <- list.files(path = out_f, pattern="RAP_TCC1km_", 
                       all.files=TRUE, full.names=FALSE)
rasters_1km <-lapply(file.path(out_f, rastlist), rast)
rasters5k<-lapply(rasters_1km, function(x) {focal(x,w=matrix(1,5,5),fun=mean,na.rm=TRUE)})
for (i in 1:length(rastlist)) {names(rasters5k[[i]]) <- paste0("RAP_TCC5x5_", substr(rastlist[i],12,15), ".tif")}
lapply(rasters5k, function(i) {
  writeRaster(i, filename=file.path(out_f,names(i)), filetype="GTiff", overwrite = TRUE)
})




  img_1k <- ee_as_raster(
     img,
      via = "drive",
      dsn = "ORNLbiomass.tif",
      #scale = 1000,
      container = "rgee_backup",
      maxPixels = 4e+09,
      lazy = FALSE,
      public = TRUE,
      add_metadata = TRUE,
      timePrefix = TRUE,
      quiet = FALSE,
   )
  return(img_1k)
}

# Run function to extract
carbonDen.extract <- carbonDen(yr)
carbonDen1k <- rast(carbonDen.extract)
carbonDen_pj <- project(carbonDen1k, rast1k, res = 1000, method = "near", align = TRUE)
carbonDen_cp <-crop(carbonDen_pj, ext(rast1k), extend = TRUE)
writeRaster(carbonDen_cp, filename=file.path(out_f, "ORNLbiomass1km.tif"), overwrite=TRUE)

writeRaster(carbonDen_pj, filename=file.path(out_f, "ORNLbiomass1km_notcrop.tif"), overwrite=TRUE)
writeRaster(carbonDen1k, filename=file.path(out_f, "ORNLbiomass1km_rast.tif"), overwrite=TRUE)



